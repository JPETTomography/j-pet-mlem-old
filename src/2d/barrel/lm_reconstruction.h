#pragma once

#include "2d/barrel/detector_set.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_map.h"
#include "2d/barrel/geometry_soa.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#define FULL_EVENT_INFO 1

namespace PET2D {
namespace Barrel {

/// 2D barrel list-mode reconstuction
template <typename FType, typename SType, std::size_t MaxDetectorsSize = 192>
class LMReconstruction {
  using Detector = PET2D::Barrel::SquareDetector<FType>;
  using Scanner2D =
      PET2D::Barrel::DetectorSet<Detector, SType, MaxDetectorsSize>;

 public:
  using F = FType;
  using S = SType;
  using Response = typename Scanner2D::Response;
  using Point = PET2D::Point<F>;
  using LOR = PET2D::Barrel::LOR<S>;
  using Grid = PET2D::PixelGrid<F, S>;
  using Geometry = PET2D::Barrel::GeometrySOA<F, S>;

  struct Event {
    LOR lor;
    F t;
#if FULL_EVENT_INFO
    Point p;
    F gauss_norm;
    F inv_sigma2;
#endif
    size_t pixel_info_begin;
    size_t pixel_info_end;
  };

#if !__CUDACC__
  using Pixel = typename Geometry::Pixel;
  using RawOutput = std::vector<F>;
  using Output = PET2D::PixelMap<Pixel, F>;

  struct PixelKernelInfo {
    S ix, iy;
    F weight;
  };

  LMReconstruction(const Grid& grid, const Geometry& geometry, F sigma)
      : grid(grid),
        geometry(geometry),
        n_pixels(grid.n_pixels),
        system_matrix_(false),
        sigma_(sigma),
        gauss_norm_dl_(1 / (sigma * std::sqrt(2 * M_PI))),
        inv_sigma2_dl_(1 / (2 * sigma * sigma)),
        rho_(grid.n_rows, grid.n_columns, 1),
        sensitivity_(grid.n_rows, grid.n_columns),
        n_threads_(omp_get_max_threads()),
        thread_rhos_(n_threads_),
        thread_kernel_caches_(n_threads_),
        n_events_per_thread_(n_threads_, 0) {}

  const Output& rho() const { return rho_; }

  void use_system_matrix() { system_matrix_ = true; }

  F sigma_w(F width) const { return 0.3 * width; }

  Event to_event(const Response& response) {
    Event event;
    event.lor = response.lor;
    const auto lor_index = event.lor.index();
    const auto& segment = geometry.lor_line_segments[lor_index];

#if FULL_EVENT_INFO
    const auto width = geometry.lor_widths[lor_index];
    event.gauss_norm = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
    event.inv_sigma2 = 1 / (2 * sigma_w(width) * sigma_w(width));
#endif
    F t = 0.5 - response.dl / (2 * segment.length);
    event.t = t;
#if FULL_EVENT_INFO
    event.p = segment.start.interpolate(segment.end, t);
#endif

    const auto t_up = t + 3 * sigma_;
    const auto t_dn = t - 3 * sigma_;

    const auto lor_info_begin = geometry.lor_pixel_info_begin[lor_index];
    const auto lor_info_end = geometry.lor_pixel_info_end[lor_index];
    const auto pixel_positions_begin =
        &geometry.pixel_positions[lor_info_begin];
    const auto pixel_positions_end = &geometry.pixel_positions[lor_info_end];

    event.pixel_info_end =
        std::upper_bound(pixel_positions_begin,
                         pixel_positions_end,
                         t_up,
                         [](const F a, const F b) -> bool { return a < b; }) -
        pixel_positions_begin + 1;

    event.pixel_info_begin =
        std::lower_bound(pixel_positions_begin,
                         pixel_positions_end,
                         t_dn,
                         [](const F a, const F b) -> bool { return a < b; }) -
        pixel_positions_begin;

#if SYSTEM_MATRIX
    if (system_matrix_) {
      event.pixel_info_end = geometry[event.lor].pixels.end();
      event.pixel_info_begin = geometry[event.lor].pixels.begin();
    }
#endif
    return event;
  }

  void add(const Response& response) { events_.push_back(to_event(response)); }

  /// Load response from input stream
  LMReconstruction& operator<<(std::istream& in) {
    for (;;) {
      Response response(in);
      if (!in)
        break;
      add(response);
    }
    return *this;
  }

  F kernel_l(const Event& event, const F pixel_position) const {
    auto segment = geometry.lor_line_segments[event.lor.index()];
    auto diff_t = (pixel_position - event.t) * segment.length;
    return gauss_norm_dl_ * compat::exp(-diff_t * diff_t * inv_sigma2_dl_);
  }

  F kernel(const Event& event,
           const Pixel pixel,
           const F pixel_position,
           const F pixel_weight) const {
    const auto pixel_index = grid.index(pixel);
    const auto kernel_z = pixel_weight / sensitivity_[pixel_index];
    return kernel_l(event, pixel_position) * kernel_z * rho_[pixel_index];
  }

  int operator()() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;

    reset_thread_data();

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // --- event loop ----------------------------------------------------------
    for (int event_index = 0; event_index < (int)events_.size();
         ++event_index) {
      int thread = omp_get_thread_num();
      n_events_per_thread_[thread]++;
      auto event = events_[event_index];

      // -- voxel loop - denominator -------------------------------------------
      F denominator = 0;
      // -- voxel loop - denominator -------------------------------------------
      for (auto info_index = event.pixel_info_begin;
           info_index < event.pixel_info_end;
           ++info_index) {
        const auto pixel = geometry.pixels[info_index];
        const auto pixel_position = geometry.pixel_positions[info_index];
        const auto pixel_weight = geometry.pixel_weights[info_index];
        const auto pixel_index = grid.index(pixel);
        const auto weight = kernel(event, pixel, pixel_position, pixel_weight);
        thread_kernel_caches_[thread][pixel_index] = weight;
        denominator += weight * sensitivity_[pixel_index];
      }  // voxel loop - denominator

      if (denominator == 0)
        continue;

      const auto inv_denominator = 1 / denominator;

      // -- voxel loop ---------------------------------------------------------
      for (auto info_index = event.pixel_info_begin;
           info_index < event.pixel_info_end;
           ++info_index) {
        const auto pixel = geometry.pixels[info_index];
        const auto pixel_index = grid.index(pixel);
        thread_rhos_[thread][pixel_index] +=
            thread_kernel_caches_[thread][pixel_index] * inv_denominator;
      }  // voxel loop
    }    // event loop
    event_count_ = 0;

    rho_.assign(0);
    for (int thread = 0; thread < n_threads_; ++thread) {
      for (int i = 0; i < grid.n_pixels; ++i) {
        rho_[i] += thread_rhos_[thread][i];
      }
      event_count_ += n_events_per_thread_[thread];
    }

    return event_count_;
  }

  void calculate_weight() {
    sensitivity_.assign(0);
    for (size_t lor_index = 0; lor_index < geometry.n_lors; ++lor_index) {
      const auto& segment = geometry.lor_line_segments[lor_index];
      const auto width = geometry.lor_widths[lor_index];
      const auto gauss_norm_w = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
      const auto inv_sigma2_w = 1 / (2 * sigma_w(width) * sigma_w(width));

      for (size_t pixel_info = geometry.lor_pixel_info_begin[lor_index];
           pixel_info < geometry.lor_pixel_info_end[lor_index];
           ++pixel_info) {
        auto pixel = geometry.pixels[pixel_info];
        auto center = grid.center_at(pixel);
        auto distance = segment.distance_from(center);
        auto kernel_z =
            gauss_norm_w * std::exp(-distance * distance * inv_sigma2_w);
        geometry.pixel_weights[pixel_info] = kernel_z;
      }
    }
  }

  void calculate_sensitivity() {
    sensitivity_.assign(0);
    for (size_t pixel_info = 0; pixel_info < geometry.n_pixel_infos;
         ++pixel_info) {
      const auto pixel = geometry.pixels[pixel_info];
      sensitivity_[pixel] += geometry.pixel_weights[pixel_info];
    }
  }

  const Output& sensitivity() const { return sensitivity_; }

  Event event(int i) const { return events_[i]; }
  const std::vector<Event>& events() const { return events_; }
  size_t n_events() const { return events_.size(); }

  F sigma() const { return sigma_; }

  const Grid grid;
  const Geometry& geometry;
  const int n_pixels;

 private:
  void reset_thread_data() {
    const auto n_pixels = grid.n_pixels;
    for (auto& thread_rho : thread_rhos_) {
      thread_rho.assign(n_pixels, 0);
    }
    for (auto& thread_kernel_cache : thread_kernel_caches_) {
      thread_kernel_cache.assign(n_pixels, 0);
    }
    for (auto& n_events : n_events_per_thread_) {
      n_events = 0;
    }
  }

  bool system_matrix_;
  std::vector<Event> events_;
  F sigma_;
  F gauss_norm_dl_;
  F inv_sigma2_dl_;

  Output rho_;
  Output sensitivity_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
  int n_threads_;

  std::vector<RawOutput> thread_rhos_;
  std::vector<RawOutput> thread_kernel_caches_;
  std::vector<int> n_events_per_thread_;
#endif  // !__CUDACC__
};

}  // Barrel
}  // PET2D

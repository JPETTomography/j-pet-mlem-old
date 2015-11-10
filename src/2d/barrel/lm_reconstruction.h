#pragma once

#include "2d/barrel/detector_set.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_map.h"
#if !__CUDACC__
#include "2d/barrel/geometry.h"
#endif

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

  struct Event {
    LOR lor;
    F t;
#if FULL_EVENT_INFO
    Point p;
    F gauss_norm;
    F inv_sigma2;
#endif
    size_t first_pixel_info_index;
    size_t last_pixel_info_index;
  };

#if !__CUDACC__
  using Geometry = PET2D::Barrel::Geometry<F, S>;
  using LORGeometry = typename Geometry::LORGeometry;
  using Pixel = typename Geometry::Pixel;
  using PixelInfo = typename Geometry::PixelInfo;
  using RawOutput = std::vector<F>;
  using Output = PET2D::PixelMap<Pixel, F>;

  struct PixelKernelInfo {
    S ix, iy;
    F weight;
  };

  LMReconstruction(Geometry& geometry, F sigma)
      : geometry(geometry),
        n_pixels(geometry.grid.n_pixels),
        system_matrix_(false),
        sigma_(sigma),
        gauss_norm_dl_(1 / (sigma * std::sqrt(2 * M_PI))),
        inv_sigma2_dl_(1 / (2 * sigma * sigma)),
        rho_(geometry.grid.n_rows, geometry.grid.n_columns, 1),
        sensitivity_(geometry.grid.n_rows, geometry.grid.n_columns),
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

    auto& segment = geometry[response.lor].segment;

#if FULL_EVENT_INFO
    auto width = geometry[event.lor].width;
    event.gauss_norm = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
    event.inv_sigma2 = 1 / (2 * sigma_w(width) * sigma_w(width));
#endif
    F t = 0.5 - response.dl / (2 * segment.length);
    event.t = t;
#if FULL_EVENT_INFO
    event.p = segment.start.interpolate(segment.end, t);
#endif

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t + 3 * sigma_;
    pix_info_dn.t = t - 3 * sigma_;
    const auto& lor_geometry = geometry[event.lor];
    event.last_pixel_info_index =
        std::upper_bound(lor_geometry.pixel_infos.begin(),
                         lor_geometry.pixel_infos.end(),
                         pix_info_up,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         }) -
        lor_geometry.pixel_infos.begin();

    event.first_pixel_info_index =
        std::lower_bound(lor_geometry.pixel_infos.begin(),
                         lor_geometry.pixel_infos.end(),
                         pix_info_dn,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         }) -
        lor_geometry.pixel_infos.begin();

#if SYSTEM_MATRIX
    if (system_matrix_) {
      event.last_pixel = geometry[event.lor].pixels.end();
      event.first_pixel = geometry[event.lor].pixels.begin();
    }
#endif
    return event;
  }

  void add(const Response& response) {
    auto event = to_event(response);
    events_.push_back(event);
  }

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

  F kernel_l(const Event& event, const PixelInfo& pixel_info) const {

    auto lor = event.lor;
    auto segment = geometry[lor].segment;
    auto diff_t = (pixel_info.t - event.t) * segment.length;
    return gauss_norm_dl_ * compat::exp(-diff_t * diff_t * inv_sigma2_dl_);
  }

  F kernel(const Event& event, const PixelInfo& pixel_info) const {
    const auto pixel = pixel_info.pixel;
    const auto pixel_index = geometry.grid.index(pixel);
    const auto kernel_z = pixel_info.weight / sensitivity_[pixel_index];
    return kernel_l(event, pixel_info) * kernel_z * rho_[pixel_index];
  }

  int operator()() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = geometry.grid;

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
      const auto& lor_geometry = geometry[event.lor];
      for (auto info_index = event.first_pixel_info_index;
           info_index < event.last_pixel_info_index;
           ++info_index) {
        pixel_count_++;
        const auto& pixel_info = lor_geometry.pixel_infos[info_index];
        const auto weight = kernel(event, pixel_info);
        thread_kernel_caches_[thread][grid.index(pixel_info.pixel)] = weight;
        denominator += weight * sensitivity_[grid.index(pixel_info.pixel)];

      }  // voxel loop - denominator

      if (denominator == 0)
        continue;

      const auto inv_denominator = 1 / denominator;

      // -- voxel loop ---------------------------------------------------------
      for (auto info_index = event.first_pixel_info_index;
           info_index < event.last_pixel_info_index;
           ++info_index) {
        const auto& pixel_info = lor_geometry.pixel_infos[info_index];
        const auto pixel_index = grid.index(pixel_info.pixel);
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
    auto& grid = geometry.grid;

    for (auto& lor_geometry : geometry) {

      auto& segment = lor_geometry.segment;

      auto width = lor_geometry.width;
      auto gauss_norm_w = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
      auto inv_sigma2_w = 1 / (2 * sigma_w(width) * sigma_w(width));

      for (auto& pixel_info : lor_geometry.pixel_infos) {
        auto pixel = pixel_info.pixel;
        auto center = grid.center_at(pixel);
        auto distance = segment.distance_from(center);
        auto kernel_z =
            gauss_norm_w * std::exp(-distance * distance * inv_sigma2_w);
        pixel_info.weight = kernel_z;
      }
    }
  }

  void calculate_sensitivity() {
    auto& grid = geometry.grid;

    sensitivity_.assign(0);
    for (auto& lor_geometry : geometry) {
      for (auto& pixel_info : lor_geometry.pixel_infos) {
        auto pixel = pixel_info.pixel;
        auto pixel_index = grid.index(pixel);
        sensitivity_[pixel_index] += pixel_info.weight;
      }
    }
  }

  const Output& sensitivity() const { return sensitivity_; }

  Event event(int i) const { return events_[i]; }
  const std::vector<Event>& events() const { return events_; }
  size_t n_events() const { return events_.size(); }

  F sigma() const { return sigma_; }

  Geometry& geometry;
  const int n_pixels;

 private:
  void reset_thread_data() {
    const auto n_pixels = geometry.grid.n_pixels;
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

#pragma once

#include "2d/barrel/geometry.h"
#include "2d/barrel/detector_set.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"
#include "2d/geometry/pixel_map.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET2D {
namespace Barrel {

/// 2D barrel list-mode reconstuction
template <typename FType, typename SType> class LMReconstruction {
  using Detector = PET2D::Barrel::SquareDetector<FType>;
  using Scanner2D = PET2D::Barrel::DetectorSet<Detector, SType>;

 public:
  using F = FType;
  using S = SType;
  using Response = typename Scanner2D::Response;
  using Point = PET2D::Point<F>;
  using LOR = PET2D::Barrel::LOR<S>;
  using Geometry = PET2D::Barrel::Geometry<F, S>;
  using LORInfo = typename Geometry::LORInfo;
  using Pixel = typename Geometry::Pixel;
  using PixelInfo = typename Geometry::PixelInfo;
  using PixelConstIterator = typename LORInfo::PixelInfoList::const_iterator;
  using RawOutput = std::vector<F>;
  using Output = PET2D::PixelMap<Pixel, F>;

  struct BarrelEvent {
    LOR lor;
    Point p;
    PixelConstIterator first_pixel;
    PixelConstIterator last_pixel;
    F t;
    F gauss_norm;
    F inv_sigma2;
  };

  struct PixelKernelInfo {
    S ix, iy;
    F weight;
  };

  LMReconstruction(Geometry& geometry, F sigma)
      : geometry(geometry),
        n_pixels(geometry.grid.n_pixels),
        system_matrix_(false),
        sigma_(sigma),
        rho_(geometry.grid.n_rows, geometry.grid.n_columns, 1),
        sensitivity_(geometry.grid.n_rows, geometry.grid.n_columns),
        n_threads_(omp_get_max_threads()),
        thread_rhos_(n_threads_),
        thread_kernel_caches_(n_threads_),
        n_events_per_thread_(n_threads_, 0) {}

  const Output& rho() const { return rho_; }

  void use_system_matrix() { system_matrix_ = true; }

  F sigma_w(F width) const { return 0.3 * width; }

  BarrelEvent to_event(const Response& response) {
    BarrelEvent event;
    event.lor = response.lor;

    auto& segment = geometry[response.lor].segment;

    auto width = geometry[event.lor].width;
    event.gauss_norm = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
    event.inv_sigma2 = 1 / (2 * sigma_w(width) * sigma_w(width));

    F t = 0.5 - response.dl / (2 * segment.length);
    event.t = t;
    event.p = segment.start.interpolate(segment.end, t);

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t + 3 * sigma_;
    pix_info_dn.t = t - 3 * sigma_;
    event.last_pixel =
        std::upper_bound(geometry[event.lor].pixels.begin(),
                         geometry[event.lor].pixels.end(),
                         pix_info_up,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         });

    event.first_pixel =
        std::lower_bound(geometry[event.lor].pixels.begin(),
                         geometry[event.lor].pixels.end(),
                         pix_info_dn,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         });

#if SYSTEM_MATRIX
    if (system_matrix_) {
      event.last_pixel = geometry[event.lor].pixels.end();
      event.first_pixel = geometry[event.lor].pixels.begin();
    }
#endif
    return event;
  }

  /// Load response from input stream
  LMReconstruction& operator<<(std::istream& in) {
    for (;;) {
      Response response(in);
      if (!in)
        break;
      auto event = to_event(response);
      events_.push_back(event);
    }
    return *this;
  }

  int operator()() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = geometry.grid;
    for (auto& thread_rho : thread_rhos_) {
      thread_rho.assign(grid.n_pixels, 0);
    }
    for (auto& thread_kernel_cache : thread_kernel_caches_) {
      thread_kernel_cache.assign(grid.n_pixels, 0);
    }
    for (auto& n_events : n_events_per_thread_) {
      n_events = 0;
    }

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // --- event loop ----------------------------------------------------------
    for (size_t i = 0; i < events_.size(); ++i) {
      int thread = omp_get_thread_num();
      n_events_per_thread_[thread]++;
      auto event = events_[i];

      // -- voxel loop - denominator -------------------------------------------
      double denominator = 0;
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        pixel_count_++;
        auto pix = it->pixel;

        int index = grid.index(pix.x, pix.y);
        double kernel_z = it->weight / sensitivity_[index];

        auto gauss_norm_dl = 1 / (sigma_ * std::sqrt(2 * M_PI));
        auto inv_sigma2_dl = 1 / (2 * sigma_ * sigma_);

        F kernel_l = gauss_norm_dl * exp(-(it->t - event.t) *
                                         (it->t - event.t) * inv_sigma2_dl);
        F weight = kernel_l * kernel_z * rho_[index];

        thread_kernel_caches_[thread][index] = weight;

        denominator += weight;

      }  // voxel loop - denominator

      double inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        continue;
      }

      // -- voxel loop ---------------------------------------------------------
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        auto pix = it->pixel;

        int index = grid.index(pix.x, pix.y);

        thread_rhos_[thread][index] +=
            thread_kernel_caches_[thread][index] * inv_denominator;

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

    for (auto& lor_info : geometry) {

      auto& segment = lor_info.segment;

      auto width = lor_info.width;
      auto gauss_norm_w = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
      auto inv_sigma2_w = 1 / (2 * sigma_w(width) * sigma_w(width));

      for (auto& pixel_info : lor_info.pixels) {
        auto pixel = pixel_info.pixel;
        auto center = grid.center_at(pixel.x, pixel.y);
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
    for (auto& lor_info : geometry) {
      for (auto& pixel_info : lor_info.pixels) {
        auto pixel = pixel_info.pixel;
        auto index = grid.index(pixel);
        sensitivity_[index] += pixel_info.weight;
      }
    }
  }

  const Output& sensitivity() const { return sensitivity_; }

  BarrelEvent event(int i) const { return events_[i]; }

  size_t n_events() const { return events_.size(); }

  Geometry& geometry;
  const int n_pixels;

 private:
  bool system_matrix_;
  std::vector<BarrelEvent> events_;
  F sigma_;

  Output rho_;
  Output sensitivity_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
  int n_threads_;

  std::vector<RawOutput> thread_rhos_;
  std::vector<RawOutput> thread_kernel_caches_;
  std::vector<int> n_events_per_thread_;
};

}  // Barrel
}  // PET2D

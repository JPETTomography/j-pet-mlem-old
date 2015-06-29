#ifndef LM_RECONSTRUCTION
#define LM_RECONSTRUCTION

#include "2d/barrel/lors_pixels_info.h"
#include "2d/barrel/detector_set.h"
#include "2d/geometry/point.h"
#include "2d/barrel/lor.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class LMReconstruction {
  using Detector = PET2D::Barrel::SquareDetector<FType>;
  using Scanner2D = PET2D::Barrel::DetectorSet<Detector, 192, SType>;

 public:
  using F = FType;
  using S = SType;
  using Response = typename Scanner2D::Response;
  using Point = PET2D::Point<F>;
  using LOR = PET2D::Barrel::LOR<S>;
  using LORsPixelsInfo = PET2D::Barrel::LORsPixelsInfo<F, S>;
  using PixelInfo = typename LORsPixelsInfo::PixelInfo;
  using PixelConstIterator =
      typename LORsPixelsInfo::PixelInfoContainer::const_iterator;

  struct BarrelEvent {
    LOR lor;
    Point p;
    PixelConstIterator first_pixel;
    PixelConstIterator last_pixel;

    F gauss_norm;
    F inv_sigma2;
  };

  struct PixelKernelInfo {
    S ix, iy;
    F weight;
  };

  LMReconstruction(const LORsPixelsInfo& lor_pixel_info, F sigma)
      : lor_pixel_info(lor_pixel_info),
        sigma_(sigma),
        n_pixels(lor_pixel_info.grid.n_pixels),
        n_threads_(omp_get_max_threads()),
        thread_rhos_(n_threads_),
        thread_kernel_caches_(n_threads_),
        n_events_per_thread_(n_threads_, 0) {}

  F sigma(F width) const { return 0.3 * width; }

  typename std::vector<F>::const_iterator rho_begin() const {
    return rho_.begin();
  }
  typename std::vector<F>::const_iterator rho_end() const { return rho_.end(); }
  typename std::vector<F>::iterator rho_begin() { return rho_.begin(); }

  BarrelEvent to_event(const Response& response) {
    BarrelEvent event;
    event.lor = response.lor;

    auto segment = lor_pixel_info[response.lor].segment;

    auto width = lor_pixel_info[event.lor].width;
    event.gauss_norm = 1 / (sigma(width) * std::sqrt(2 * M_PI));
    event.inv_sigma2 = 1 / (2 * sigma(width) * sigma(width));

    F t = 0.5 - response.dl / (2 * segment->length);
    event.p = PET2D::interpolate(t, segment->start, segment->end);

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t + 3 * sigma_;
    pix_info_dn.t = t - 3 * sigma_;
    event.last_pixel =
        std::upper_bound(lor_pixel_info[event.lor].pixels.begin(),
                         lor_pixel_info[event.lor].pixels.end(),
                         pix_info_up,
                         [](const PixelInfo& a, const PixelInfo& b)
                             -> bool { return a.t < b.t; });

    event.first_pixel =
        std::lower_bound(lor_pixel_info[event.lor].pixels.begin(),
                         lor_pixel_info[event.lor].pixels.end(),
                         pix_info_dn,
                         [](const PixelInfo& a, const PixelInfo& b)
                             -> bool { return a.t < b.t; });
    return event;
  }

  int fscanf_responses(std::istream& in) {
    int count = 0;
    while (true) {
      auto response = fscanf_response(in);
      if (!in)
        break;
      count++;
      auto event = to_event(response);
      events_.push_back(event);
    }
    return count;
  }

  int iterate() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = lor_pixel_info.grid;
    for (auto& thread_rho : thread_rhos_) {
      thread_rho.assign(grid.n_pixels, 0);
    }
    for (auto& thread_kernel_cache : thread_kernel_caches_) {
      thread_kernel_cache.assign(grid.n_pixels, 0);
    }

    for (auto& n_events : n_events_per_thread_) {
      n_events = 0;
    }

/* ------- Event loop ------*/
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < events_.size(); ++i) {
      int thread = omp_get_thread_num();
      n_events_per_thread_[thread]++;
      auto event = events_[i];
      auto lor = event.lor;
      auto segment = *lor_pixel_info[lor].segment;
      auto R = segment.length / 2;

      /* ---------  Voxel loop  - denominator ----------- */
      double denominator = 0;
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        pixel_count_++;
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;

        auto center = grid.center_at(ix, iy);
        auto distance = segment.distance_from(center);

        int index = grid.index(ix, iy);

        auto kernel_z = event.gauss_norm *
                        std::exp(-distance * distance * event.inv_sigma2);
        auto weight = kernel_z * rho_[index];

        thread_kernel_caches_[thread][index] = weight;
        denominator += weight;

      }  // Voxel loop - denominator

      F inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        std::cerr << "denminator == 0 !";
        abort();
      }

      /* ---------  Voxel loop ------------ */
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;

        int index = grid.index(ix, iy);
        thread_rhos_[thread][index] +=
            thread_kernel_caches_[thread][index] * inv_denominator;

      }  // Voxel loop
    }    // event loop
    event_count_ = 0;

    rho_.assign(grid.n_pixels, 0);
    for (int thread = 0; thread < n_threads_; ++thread) {
      for (int i = 0; i < grid.n_pixels; ++i) {
        rho_[i] += thread_rhos_[thread][i];
      }
      event_count_ += n_events_per_thread_[thread];
    }

    return event_count_;
  }

  const LORsPixelsInfo& lor_pixel_info;
  const int n_pixels;

 private:
  Response fscanf_response(std::istream& in) {
    S d1, d2;
    Response response;
    in >> d1 >> d2 >> response.tof_position >> response.dl;
    response.lor = LOR(d1, d2);
    return response;
  }

  std::vector<BarrelEvent> events_;
  F sigma_;

  std::vector<F> rho_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
  int n_threads_;

  std::vector<std::vector<F>> thread_rhos_;
  std::vector<std::vector<F>> thread_kernel_caches_;
  std::vector<PixelKernelInfo> voxel_cache_;
  std::vector<int> n_events_per_thread_;
};
}
}
#endif  // LM_RECONSTRUCTION

#pragma once

#include <vector>
#include <algorithm>

#include "2d/barrel/geometry.h"
#include "2d/strip/response.h"
#include "3d/geometry/point.h"
#include "3d/geometry/voxel_grid.h"

#include "common/mathematica_graphics.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET3D {
namespace Hybrid {

/// 3D hybrid PET reconstruction
template <class ScannerClass, class Kernel2DClass> class Reconstruction {
 public:
  using Scanner = ScannerClass;
  using Kernel2D = Kernel2DClass;
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using Geometry = PET2D::Barrel::Geometry<F, S>;
  using LORInfo = typename Geometry::LORInfo;
  using Response = typename Scanner::Response;
  using LOR = PET2D::Barrel::LOR<S>;
  using StripEvent = PET2D::Strip::Response<F>;
  using PixelInfo = typename Geometry::PixelInfo;
  using Pixel = typename Geometry::Pixel;
  using Point2D = PET2D::Point<F>;
  using Point = PET3D::Point<F>;
  using Vector2D = PET2D::Vector<F>;
  using PixelConstIterator = typename LORInfo::PixelInfoList::const_iterator;
  using Output = std::vector<F>;

  struct FrameEvent {
    LOR lor;
    F up;
    F right;
    F tan;
    F sec;
    F half_box_up;
    F half_box_right;
    PixelConstIterator first_pixel;
    PixelConstIterator last_pixel;
    S first_plane;
    S last_plane;
    F gauss_norm;
    F inv_sigma2;
  };

  struct VoxelKernelInfo {
    S ix, iy, iz;
    F weight;
  };

  Reconstruction(const Scanner& scanner,
                 Geometry& geometry,
                 F z_left,
                 int n_planes)
      : scanner(scanner),
        geometry(geometry),
        z_left(z_left),
        n_planes(n_planes),
        v_grid(geometry.grid, z_left, n_planes),
        n_voxels(geometry.grid.n_columns * geometry.grid.n_rows * n_planes),
        kernel_(scanner.sigma_z(), scanner.sigma_dl()),
        rho_(n_voxels, F(1.0)),
        n_threads_(omp_get_max_threads()),
        thread_rhos_(n_threads_),
        thread_kernel_caches_(n_threads_),
        n_events_per_thread_(n_threads_, 0) {}

  Point translate_to_point(const Response& response) {

    auto segment = geometry[response.lor].segment;
    F t = 0.5 - response.dl / (2 * segment->length);
    return Point(segment->start.x, segment->start.y, response.z_dn)
        .iterpolate(Point(segment->end.x, segment->end.y, response.z_up), t);
  }

  F sigma_w(F width) const { return 0.3 * width; }

  FrameEvent translate_to_frame(const Response& response) {

    FrameEvent event;
    event.lor = response.lor;

    auto R = geometry[event.lor].segment.length / 2;
    StripEvent strip_event(response.z_up, response.z_dn, response.dl);

    auto width = geometry[event.lor].width;
    event.gauss_norm = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
    event.inv_sigma2 = 1 / (2 * sigma_w(width) * sigma_w(width));
    strip_event.transform(R, event.tan, event.up, event.right);
    F A, B, C;
    kernel_.ellipse_bb(
        event.tan, event.sec, A, B, C, event.half_box_up, event.half_box_right);

    auto ev_z_left = event.right - event.half_box_right;
    auto ev_z_right = event.right + event.half_box_right;
    event.first_plane = plane(ev_z_left);
    event.last_plane = plane(ev_z_right) + 1;

    auto y_up = event.up + event.half_box_up;
    auto y_dn = event.up - event.half_box_up;
    auto t_up = (y_up + R) / (2 * R);
    auto t_dn = (y_dn + R) / (2 * R);

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t_up;
    pix_info_dn.t = t_dn;
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

    return event;
  }

  S plane(F z) {
    return S(std::floor((z - z_left) / geometry.grid.pixel_size));
  }

  Reconstruction& operator<<(std::istream& in) {
    for (;;) {
      Response response(in);
      if (!in)
        break;
      events_.push_back(translate_to_frame(response));
    }
    return *this;
  }

  int n_events() const { return events_.size(); }
  FrameEvent frame_event(int i) const { return events_[i]; }

  void calculate_weight() {
    auto& grid = geometry.grid;
    sensitivity_.assign(grid.n_pixels, 0);
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
    sensitivity_.assign(grid.n_pixels, 0);
    for (auto& lor_info : geometry) {

      for (auto& pixel_info : lor_info.pixels) {
        auto pixel = pixel_info.pixel;
        auto index = grid.index(pixel);
        sensitivity_[index] += pixel_info.weight;
      }
    }
  }

  int iterate() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = geometry.grid;
    for (auto& thread_rho : thread_rhos_) {
      thread_rho.assign(v_grid.n_voxels, 0);
    }
    for (auto& thread_kernel_cache : thread_kernel_caches_) {
      thread_kernel_cache.assign(v_grid.n_voxels, 0);
    }

    for (auto& n_events : n_events_per_thread_) {
      n_events = 0;
    }

/* ------- Event loop ------*/
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n_events(); ++i) {
      int thread = omp_get_thread_num();
      n_events_per_thread_[thread]++;
      auto event = frame_event(i);
      auto lor = event.lor;
      auto& segment = geometry[lor].segment;
      auto R = segment.length / 2;

      /* ---------  Voxel loop  - denominator ----------- */
      double denominator = 0;
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        pixel_count_++;
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;

        auto center = grid.center_at(ix, iy);
        auto up = segment.projection_relative_middle(center);

        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          voxel_count_++;
          auto z = v_grid.center_z_at(ix, iy, iz);
          int index = v_grid.index(ix, iy, iz);

          auto diff = Point2D(up, z) - Point2D(event.up, event.right);
          auto kernel2d = kernel_(
              event.up, event.tan, event.sec, R, Vector2D(diff.y, diff.x));
          auto kernel_z = it->weight;
          auto weight = kernel2d * kernel_z * rho_[index] /
                        sensitivity_[grid.index(ix, iy)];

          thread_kernel_caches_[thread][index] = weight;
          denominator += weight;
        }
      }  // Voxel loop - denominator

      F inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        throw("denminator == 0 !");
      }

      /* ---------  Voxel loop ------------ */
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        auto pix = it->pixel;
        auto ix = pix.x;
        auto iy = pix.y;
        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          int index = v_grid.index(ix, iy, iz);
          thread_rhos_[thread][index] +=
              thread_kernel_caches_[thread][index] * inv_denominator;
        }
      }  // Voxel loop
    }    // event loop
    event_count_ = 0;

    rho_.assign(v_grid.n_voxels, 0);
    for (int thread = 0; thread < n_threads_; ++thread) {
      for (int i = 0; i < v_grid.n_voxels; ++i) {
        rho_[i] += thread_rhos_[thread][i];
      }
      event_count_ += n_events_per_thread_[thread];
    }
    // std::swap(rho_, rho_new_);
    return event_count_;
  }

  const Output& rho() const { return rho_; }

  int voxel_count() const { return voxel_count_; }
  int pixel_count() const { return pixel_count_; }
  int event_count() const { return event_count_; }

  void graph_frame_event(Common::MathematicaGraphics<F>& graphics,
                         int event_index) {
    auto event = events_[event_index];
    auto lor = event.lor;
    graphics.add(scanner.barrel, lor);
    graphics.add(geometry[lor].segment);
    for (auto pix = event.first_pixel; pix != event.last_pixel; ++pix) {
      graphics.add_pixel(geometry.grid, pix->pixel);
    }
  }

 public:
  const Scanner& scanner;
  Geometry& geometry;
  const F z_left;
  const S n_planes;
  const PET3D::VoxelGrid<F, S> v_grid;
  const int n_voxels;

 private:
  std::vector<FrameEvent> events_;
  Kernel2D kernel_;
  Output rho_;
  int event_count_;
  int voxel_count_;
  int pixel_count_;
  int n_threads_;
  std::vector<Output> thread_rhos_;
  std::vector<Output> thread_kernel_caches_;
  std::vector<VoxelKernelInfo> voxel_cache_;
  std::vector<int> n_events_per_thread_;
  Output sensitivity_;
};

}  // Hybrid
}  // PET3D

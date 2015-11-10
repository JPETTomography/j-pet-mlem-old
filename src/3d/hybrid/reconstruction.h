#pragma once

#if !__CUDACC__
#include <vector>
#include <algorithm>

#include "common/mathematica_graphics.h"
#include "2d/barrel/geometry.h"
#endif

#include "2d/strip/response.h"
#include "3d/geometry/point.h"
#include "3d/geometry/voxel_grid.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"

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
  using Response = typename Scanner::Response;
  using LOR = PET2D::Barrel::LOR<S>;
  using StripEvent = PET2D::Strip::Response<F>;
  using Voxel = PET3D::Voxel<S>;
  using Point2D = PET2D::Point<F>;
  using Point = PET3D::Point<F>;
  using Vector2D = PET2D::Vector<F>;
  using Output = PET3D::VoxelMap<Voxel, F>;

  struct FrameEvent {
    LOR lor;
    F up;
    F right;
    F tan;
    F sec;
    size_t first_pixel_info_index;
    size_t last_pixel_info_index;
    S first_plane;
    S last_plane;
  };

  struct VoxelKernelInfo {
    S ix, iy, iz;
    F weight;
  };

#if !__CUDACC__
  using Geometry = PET2D::Barrel::Geometry<F, S>;
  using LORGeometry = typename Geometry::LORGeometry;
  using PixelInfo = typename Geometry::PixelInfo;
  using Pixel = typename Geometry::Pixel;

  Reconstruction(const Scanner& scanner,
                 Geometry& geometry,
                 F z_left,
                 int n_planes)
      : scanner(scanner),
        geometry(geometry),
        z_left(z_left),
        n_planes(n_planes),
        v_grid(geometry.grid, z_left, n_planes),
        kernel_(scanner.sigma_z(), scanner.sigma_dl()),
        rho_(geometry.grid.n_columns, geometry.grid.n_rows, n_planes, 1),
        n_threads_(omp_get_max_threads()),
        n_events_per_thread_(n_threads_, 0),
        sensitivity_(geometry.grid.n_columns, geometry.grid.n_rows, n_planes) {
    for (int i = 0; i < n_threads_; ++i) {
      thread_rhos_.emplace_back(
          geometry.grid.n_columns, geometry.grid.n_rows, n_planes);
      thread_kernel_caches_.emplace_back(
          geometry.grid.n_columns, geometry.grid.n_rows, n_planes);
    }
  }

  F sigma_w(F width) const { return F(0.3) * width; }

  Point translate_to_point(const Response& response) {

    auto segment = geometry[response.lor].segment;
    F t = F(0.5) - response.dl / (2 * segment->length);
    return Point(segment->start.x, segment->start.y, response.z_dn)
        .iterpolate(Point(segment->end.x, segment->end.y, response.z_up), t);
  }

  FrameEvent translate_to_frame(const Response& response) {

    FrameEvent event;
    event.lor = response.lor;

    auto R = geometry[event.lor].segment.length / 2;
    StripEvent strip_event(response.z_up, response.z_dn, response.dl);

    strip_event.calculate_tan_y_z(R, event.tan, event.up, event.right);

    F A, B, C;
    F half_box_up, half_box_right;
    kernel_.ellipse_bb(
        event.tan, event.sec, A, B, C, half_box_up, half_box_right);

    auto ev_z_left = event.right - half_box_right;
    auto ev_z_right = event.right + half_box_right;
    event.first_plane = std::max((short)0, plane(ev_z_left));
    event.last_plane = std::min(plane(ev_z_right) + 1, v_grid.n_planes - 1);

    auto y_up = event.up + half_box_up;
    auto y_dn = event.up - half_box_up;
    auto t_up = (y_up + R) / (2 * R);
    auto t_dn = (y_dn + R) / (2 * R);

    PixelInfo pix_info_up, pix_info_dn;
    pix_info_up.t = t_up;
    pix_info_dn.t = t_dn;
    event.last_pixel_info_index =
        std::upper_bound(geometry[event.lor].pixel_infos.begin(),
                         geometry[event.lor].pixel_infos.end(),
                         pix_info_up,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         }) -
        geometry[event.lor].pixel_infos.begin();

    event.first_pixel_info_index =
        std::lower_bound(geometry[event.lor].pixel_infos.begin(),
                         geometry[event.lor].pixel_infos.end(),
                         pix_info_dn,
                         [](const PixelInfo& a, const PixelInfo& b) -> bool {
                           return a.t < b.t;
                         }) -
        geometry[event.lor].pixel_infos.begin();

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
  const std::vector<FrameEvent>& events() const { return events_; }
  FrameEvent frame_event(int i) const { return events_[i]; }

  void calculate_weight() {
    auto& grid = geometry.grid;
    sensitivity_.assign(0);
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
        auto index = grid.index(pixel);
        sensitivity_[index] += pixel_info.weight;
      }
    }
  }

  void normalize_geometry_weights() {
    auto& grid = geometry.grid;

    for (auto& lor_geometry : geometry) {

      for (auto& pixel_info : lor_geometry.pixel_infos) {
        auto pixel = pixel_info.pixel;
        auto index = grid.index(pixel);
        if (sensitivity_[index] > 0)
          pixel_info.weight /= sensitivity_[index];
      }
    }
  }

  int operator()() {
    event_count_ = 0;
    voxel_count_ = 0;
    pixel_count_ = 0;
    auto grid = geometry.grid;
    for (auto& thread_rho : thread_rhos_) {
      thread_rho.assign(0);
    }
    for (auto& thread_kernel_cache : thread_kernel_caches_) {
      thread_kernel_cache.assign(0);
    }
    for (auto& n_events : n_events_per_thread_) {
      n_events = 0;
    }

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    // --- event loop ----------------------------------------------------------
    for (int i = 0; i < n_events(); ++i) {
      int thread = omp_get_thread_num();
      n_events_per_thread_[thread]++;
      auto event = frame_event(i);
      auto lor = event.lor;
      auto& segment = geometry[lor].segment;
      auto R = segment.length / 2;

      // -- voxel loop - denominator -------------------------------------------
      double denominator = 0;
      const auto& lor_geometry = geometry[event.lor];
      for (auto info_index = event.first_pixel_info_index;
           info_index < event.last_pixel_info_index;
           ++info_index) {
        pixel_count_++;
        const auto& pixel_info = lor_geometry.pixel_infos[info_index];
        auto pixel = pixel_info.pixel;
        auto index = grid.index(pixel);
        auto center = grid.center_at(pixel);
        auto up = segment.projection_relative_middle(center);

        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          voxel_count_++;
          Voxel voxel(pixel.x, pixel.y, iz);
          auto z = v_grid.center_z_at(voxel);
          int v_index = v_grid.index(voxel);

          auto diff = Point2D(up, z) - Point2D(event.up, event.right);
          auto kernel2d = kernel_(
              event.up, event.tan, event.sec, R, Vector2D(diff.y, diff.x));
          auto kernel_t = pixel_info.weight;

          auto weight = kernel2d * kernel_t * rho_[v_index];
          denominator += weight * sensitivity_[index];

          thread_kernel_caches_[thread][v_index] = weight;
        }
      }  // voxel loop - denominator

      F inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        throw("denminator == 0 !");
      }

      // -- voxel loop ---------------------------------------------------------
      for (auto info_index = event.first_pixel_info_index;
           info_index < event.last_pixel_info_index;
           ++info_index) {
        const auto& pixel_info = lor_geometry.pixel_infos[info_index];
        auto pixel = pixel_info.pixel;
        for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
          Voxel voxel(pixel.x, pixel.y, iz);
          int index = v_grid.index(voxel);

          thread_rhos_[thread][index] +=
              thread_kernel_caches_[thread][index] * inv_denominator;
        }
      }  // voxel loop
    }    // event loop
    event_count_ = 0;

    rho_.assign(0);
    for (int thread = 0; thread < n_threads_; ++thread) {
      for (int i = 0; i < v_grid.n_voxels; ++i) {
        rho_[i] += thread_rhos_[thread][i];
      }
      event_count_ += n_events_per_thread_[thread];
    }
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
    for (auto info_index = event.first_pixel_info_index;
         info_index < event.last_pixel_info_index;
         ++info_index) {
      const auto& pixel_info = geometry[lor].pixel_infos[info_index];
      graphics.add_pixel(geometry.grid, pixel_info.pixel);
    }
  }

 public:
  const Scanner& scanner;
  Geometry& geometry;
  const F z_left;
  const S n_planes;
  const PET3D::VoxelGrid<F, S> v_grid;

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
#endif  // !__CUDACC__
};

}  // Hybrid
}  // PET3D

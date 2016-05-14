#pragma once

#if !__CUDACC__
#include <vector>
#include <algorithm>

#include "common/mathematica_graphics.h"
#if USE_FAST_TEXT_PARSER
#include "util/text_parser.h"
#endif
#endif

#include "2d/strip/response.h"
#include "2d/barrel/geometry_soa.h"
#include "2d/geometry/pixel_map.h"
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
  using NaiveOutput = PET3D::VoxelMap<Voxel, int>;
  using Grid = PET3D::VoxelGrid<F, S>;
  using PixelGrid = typename Grid::PixelGrid;
  using Pixel = typename PixelGrid::Pixel;
  using Map2D = PET2D::PixelMap<Pixel, F>;
  using Geometry = PET2D::Barrel::GeometrySOA<F, S>;

  struct FrameEvent {
    LOR lor;
    F up;
    F right;
    F tan;
    F sec;
    size_t pixel_info_begin;
    size_t pixel_info_end;
    S plane_begin;
    S plane_end;
  };

  struct VoxelKernelInfo {
    S ix, iy, iz;
    F weight;
  };

#if !__CUDACC__
  Reconstruction(const Scanner& scanner,
                 const Grid& grid,
                 const Geometry& geometry,
                 bool use_3d_sensitivity = true)
      : scanner(scanner),
        grid(grid),
        geometry(geometry),
        rho(grid.pixel_grid.n_columns,
            grid.pixel_grid.n_rows,
            grid.n_planes,
            1),
        sensitivity(grid.pixel_grid.n_columns,
                    grid.pixel_grid.n_rows,
                    geometry.n_planes_half > 1
                        ? geometry.n_planes_half
                        : use_3d_sensitivity ? (grid.n_planes / 2) : 1),
        kernel_(scanner.sigma_z(), scanner.sigma_dl()),
        n_threads_(omp_get_max_threads()),
        n_events_per_thread_(n_threads_, 0) {}

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
    const auto lor_index = event.lor.index();
    const auto& segment = geometry.lor_line_segments[lor_index];
    const auto R = segment.length / 2;
    StripEvent strip_event(response.z_up, response.z_dn, response.dl);

    strip_event.calculate_tan_y_z(R, event.tan, event.up, event.right);

    F A, B, C;
    F half_box_up, half_box_right;
    kernel_.ellipse_bb(
        event.tan, event.sec, A, B, C, half_box_up, half_box_right);

    auto ev_z_left = event.right - half_box_right;
    auto ev_z_right = event.right + half_box_right;
    event.plane_begin = std::max((S)0, plane(ev_z_left));
    event.plane_end = std::min((S)(plane(ev_z_right) + 1), grid.n_planes);
    auto y_up = event.up + half_box_up;
    auto y_dn = event.up - half_box_up;
    auto t_up = (y_up + R) / (2 * R);
    auto t_dn = (y_dn + R) / (2 * R);

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
        geometry.pixel_positions + 1;

    event.pixel_info_begin =
        std::lower_bound(pixel_positions_begin,
                         pixel_positions_end,
                         t_dn,
                         [](const F a, const F b) -> bool { return a < b; }) -
        geometry.pixel_positions;

    return event;
  }

  S plane(F z) { return S((z - grid.z_left) / grid.pixel_grid.pixel_size); }

  bool bb_intersects_grid_with_positive_weight(const FrameEvent& event) {
    if (event.plane_end < 0 || event.plane_begin >= grid.n_planes)
      return false;

    for (size_t i = event.pixel_info_begin; i < event.pixel_info_end; i++) {
      auto pixel = geometry.pixels[i];
      if (grid.pixel_grid.contains(pixel) && geometry.pixel_weights[i] > 0)
        return true;
    }

    return false;
  }

#if USE_FAST_TEXT_PARSER
  void fast_load_txt_events(const char* fn) {
    size_t n_lines = 0;
    // first just count lines and reserve space
    util::text_parser::read_lines(fn, [&](const char*) { ++n_lines; });
    events_.reserve(n_lines);
    // now read actual values
    util::text_parser::read_lines(
        fn,
        [&](const char* line) {
          util::text_parser parser(line);
          Response response;
          try {
            parser >> response.lor.first >> response.lor.second >>
                response.z_up >> response.z_dn >> response.dl;
          } catch (const char* ex) {
            std::cerr << "error line: " << line << std::endl;
            throw(ex);
          }
          auto event = translate_to_frame(response);
          if (bb_intersects_grid_with_positive_weight(event))
            events_.push_back(event);
        });
  }
#endif

  Reconstruction& operator<<(std::istream& in) {
    for (;;) {
      Response response(in);
      if (!in)
        break;
      auto event = translate_to_frame(response);
      if (bb_intersects_grid_with_positive_weight(event))
        events_.push_back(event);
    }
    return *this;
  }

  Reconstruction& operator<<(util::ibstream& in) {
    for (;;) {
      Response response(in);
      if (!in)
        break;
      auto event = translate_to_frame(response);
      if (bb_intersects_grid_with_positive_weight(event))
        events_.push_back(event);
    }
    return *this;
  }

  int n_events() const { return events_.size(); }
  const std::vector<FrameEvent>& events() const { return events_; }
  FrameEvent frame_event(int i) const { return events_[i]; }

  void calculate_weight() {
    const auto& pixel_grid = grid.pixel_grid;
    sensitivity.assign(0);
    for (size_t lor_index = 0; lor_index < geometry.n_lors; ++lor_index) {
      const auto& segment = geometry.lor_line_segments[lor_index];
      const auto width = geometry.lor_widths[lor_index];
      const auto gauss_norm_w = 1 / (sigma_w(width) * std::sqrt(2 * M_PI));
      const auto inv_sigma2_w = 1 / (2 * sigma_w(width) * sigma_w(width));

      for (size_t pixel_info = geometry.lor_pixel_info_begin[lor_index];
           pixel_info < geometry.lor_pixel_info_end[lor_index];
           ++pixel_info) {
        auto pixel = geometry.pixels[pixel_info];
        auto center = pixel_grid.center_at(pixel);
        auto distance = segment.distance_from(center);
        auto kernel_z =
            gauss_norm_w * std::exp(-distance * distance * inv_sigma2_w);
        geometry.pixel_weights[pixel_info] = kernel_z;
      }
    }
  }

  void calculate_sensitivity() {
    sensitivity.assign(0);
    for (size_t pixel_info = 0; pixel_info < geometry.n_pixel_infos;
         ++pixel_info) {
      const auto pixel = geometry.pixels[pixel_info];
      for (size_t plane = 0; plane < geometry.n_planes_half; ++plane) {
        Voxel voxel(pixel.x, pixel.y, plane);
        const auto voxel_index = pixel_info + geometry.n_pixel_infos * plane;
        sensitivity[voxel] += geometry.pixel_weights[voxel_index];
      }
    }
  }

  void normalize_geometry_weights() {
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (size_t pixel_info = 0; pixel_info < geometry.n_pixel_infos;
         ++pixel_info) {
      const auto pixel = geometry.pixels[pixel_info];
      for (size_t plane = 0; plane < geometry.n_planes_half; ++plane) {
        Voxel voxel(pixel.x, pixel.y, plane);
        const auto voxel_index = pixel_info + geometry.n_pixel_infos * plane;
        geometry.pixel_weights[voxel_index] /= sensitivity[voxel];
      }
    }
  }

  void set_sensitivity_to_one() { sensitivity.assign(1); }

  int operator()() {
    bool multiplane = geometry.n_planes_half > 1;
    bool use_3d_sensitivity = sensitivity.depth > 1;

    if (thread_rhos_.size() == 0) {
      for (int i = 0; i < n_threads_; ++i) {
        thread_rhos_.emplace_back(
            grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes);
        thread_kernel_caches_.emplace_back(
            grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes);
      }
    }

    size_t used_pixels = 0, used_voxels = 0, used_events = 0;

    const auto& pixel_grid = grid.pixel_grid;
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
      const auto event = frame_event(i);
      const auto lor = event.lor;
      const auto lor_index = lor.index();
      const auto& segment = geometry.lor_line_segments[lor_index];
      const auto R = segment.length / 2;
      F denominator = 0;

      // -- voxel loop - denominator -------------------------------------------
      for (auto info_index = event.pixel_info_begin;
           info_index < event.pixel_info_end;
           ++info_index) {
        used_pixels++;  // statistics
        const auto pixel = geometry.pixels[info_index];
        const auto pixel_weight = geometry.pixel_weights[info_index];

        const auto pixel_index = pixel_grid.index(pixel);
        const auto center = pixel_grid.center_at(pixel);
        const auto up = segment.projection_relative_middle(center);

        for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
          used_voxels++;  // statistics
          const Voxel voxel(pixel.x, pixel.y, iz);
          const auto z = grid.center_z_at(voxel);
          const auto voxel_index = grid.index(voxel);
          const auto kernel2d =
              kernel_.normalized(Point2D(event.right, event.up),
                                 event.tan,
                                 event.sec,
                                 R,
                                 scanner.length,
                                 Point2D(z, up));
          if (multiplane) {
            const auto abs_plane =
                compat::abs(iz - (int)geometry.n_planes_half);
            const auto kernel_t =
                geometry.pixel_weights[abs_plane * geometry.n_pixel_infos +
                                       info_index];
            const auto weight = kernel2d * kernel_t * rho[voxel_index];
            const auto abs_voxel = Voxel(pixel.x, pixel.y, abs_plane);
            denominator += weight * sensitivity[abs_voxel];
            thread_kernel_caches_[thread][voxel_index] = weight;
          } else if (use_3d_sensitivity) {
            const auto abs_plane = compat::abs(iz - (int)sensitivity.depth);
            const auto kernel_t = pixel_weight;
            const auto weight = kernel2d * kernel_t * rho[voxel_index];
            const auto abs_voxel = Voxel(pixel.x, pixel.y, abs_plane);
            denominator += weight * sensitivity[abs_voxel];
            thread_kernel_caches_[thread][voxel_index] = weight;
          } else {
            const auto kernel_t = pixel_weight;
            const auto weight = kernel2d * kernel_t * rho[voxel_index];
            denominator += weight * sensitivity[pixel_index];
            thread_kernel_caches_[thread][voxel_index] = weight;
          }
        }
      }  // voxel loop - denominator

      F inv_denominator;
      if (denominator > 0) {
        inv_denominator = 1 / denominator;
      } else {
        std::cerr << "denominator == < 0 !\n";
        std::cerr << "event " << i << "\n";
        std::cerr << "denominator " << denominator << "\n";
        std::cerr << "planes " << event.plane_begin << ":" << event.plane_end
                  << "\n";
        std::cerr << "(" << geometry.pixels[event.pixel_info_begin].x << ","
                  << geometry.pixels[event.pixel_info_begin].y << ") to ("
                  << geometry.pixels[event.pixel_info_end].x << ","
                  << geometry.pixels[event.pixel_info_end].y << ")\n";
        std::cerr << "emission " << point(event) << "\n";
        std::cerr << "lor " << event.lor.first << " " << event.lor.second
                  << "\n";

        for (auto info_index = event.pixel_info_begin;
             info_index < event.pixel_info_end;
             ++info_index) {

          const auto pixel = geometry.pixels[info_index];
          const auto pixel_weight = geometry.pixel_weights[info_index];
          const auto center = pixel_grid.center_at(pixel);

          std::cerr << "(" << pixel.x << "," << pixel.y << ") " << pixel_weight
                    << " " << center << "\n";
        }

        throw("denominator == 0 !");
      }

      // -- voxel loop ---------------------------------------------------------
      for (auto info_index = event.pixel_info_begin;
           info_index < event.pixel_info_end;
           ++info_index) {
        const auto pixel = geometry.pixels[info_index];
        for (auto iz = event.plane_begin; iz < event.plane_end; ++iz) {
          const Voxel voxel(pixel.x, pixel.y, iz);
          const auto voxel_index = grid.index(voxel);

          thread_rhos_[thread][voxel_index] +=
              thread_kernel_caches_[thread][voxel_index] * inv_denominator;
        }
      }  // voxel loop
    }    // event loop

    rho.assign(0);
    for (int thread = 0; thread < n_threads_; ++thread) {
      for (int i = 0; i < grid.n_voxels; ++i) {
        rho[i] += thread_rhos_[thread][i];
      }
      used_events += n_events_per_thread_[thread];
    }

    // save statistics
    statistics_.used_pixels = used_pixels;
    statistics_.used_voxels = used_voxels;
    statistics_.used_events = used_events;

    return used_events;
  }

  Point point(const FrameEvent& event) {
    const auto& segment = geometry.lor_line_segments[event.lor.index()];
    const auto point2d = segment.mid_point + segment.direction * event.up;
    return Point(point2d.x, point2d.y, event.right);
  }

  NaiveOutput naive() {
    NaiveOutput image(
        grid.pixel_grid.n_columns, grid.pixel_grid.n_rows, grid.n_planes, 0);
    for (const auto& event : events_) {
      auto p = point(event);
      const auto voxel = grid.voxel_at(p);
      if (grid.contains(voxel)) {
        ++image[voxel];
      }
    }
    return image;
  }

  void graph_frame_event(Common::MathematicaGraphics<F>& graphics,
                         int event_index) {
    auto event = events_[event_index];
    auto lor = event.lor;
    graphics.add(scanner.barrel, lor);
    graphics.add(geometry[lor].segment);
    for (auto info_index = event.pixel_info_begin;
         info_index < event.pixel_info_end;
         ++info_index) {
      const auto& pixel_info = geometry[lor].pixel_infos[info_index];
      graphics.add_pixel(grid.pixel_grid, pixel_info.pixel);
    }
  }

  /// Event statistics
  struct EventStatistics {
    size_t min_pixels, max_pixels;  ///< min/max of 2D barrel pixel number
    size_t min_planes, max_planes;  ///< min/max of 2D strip plane number
    size_t min_voxels, max_voxels;  ///< min/max of voxel number
    double avg_pixels, avg_planes, avg_voxels;  ///< averages
  };

  /// Calculates event statistics
  void event_statistics(EventStatistics& st) {
    st.min_pixels = st.min_planes = st.min_voxels = grid.n_voxels;
    st.max_pixels = st.max_planes = st.max_voxels = 0;
    size_t total_pixels = 0, total_planes = 0, total_voxels = 0;
    for (const auto& event : events_) {
      auto pixels = event.pixel_info_end - event.pixel_info_begin;
      size_t planes = event.plane_end - event.plane_begin;
      auto voxels = pixels * planes;
      total_pixels += pixels;
      total_planes += planes;
      total_voxels += voxels;
      if (pixels < st.min_pixels)
        st.min_pixels = pixels;
      if (planes < st.min_planes)
        st.min_planes = planes;
      if (voxels < st.min_voxels)
        st.min_voxels = voxels;
      if (pixels > st.max_pixels)
        st.max_pixels = pixels;
      if (planes > st.max_planes)
        st.max_planes = planes;
      if (voxels > st.max_voxels)
        st.max_voxels = voxels;
    }
    st.avg_pixels = (double)total_pixels / n_events();
    st.avg_planes = (double)total_planes / n_events();
    st.avg_voxels = (double)total_voxels / n_events();
  }

  /// Reconstruction statistics
  struct Statistics {
    size_t used_pixels;  ///< number of pixels used for reconstruction
    size_t used_voxels;  ///< number of voxels used for reconstruction
    size_t used_events;  ///< number of events used for reconstruction
  };

  /// Return reconstruction statistics
  const Statistics& statistics() const { return statistics_; }

 public:
  const Scanner& scanner;
  const Grid grid;
  const Geometry& geometry;
  Output rho;
  Output sensitivity;

 private:
  std::vector<FrameEvent> events_;
  Kernel2D kernel_;
  Statistics statistics_;
  int n_threads_;
  std::vector<Output> thread_rhos_;
  std::vector<Output> thread_kernel_caches_;
  std::vector<VoxelKernelInfo> voxel_cache_;
  std::vector<int> n_events_per_thread_;
#endif  // !__CUDACC__
};

}  // Hybrid
}  // PET3D

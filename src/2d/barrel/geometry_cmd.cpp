/// \page cmd_2d_barrel_geometry 2d_barrel_geometry
/// \brief 2D Barrel PET geometry description construction tool
///
/// Creates geometry description binary file for LM reconstruction
/// \ref cmd_2d_barrel_lm_reconstruction.
///
/// This is alternative for \ref cmd_2d_barrel_matrix. It does not use
/// Monte-Carlo, but calculates every LOR geometry and pixels that lie inside
/// this LOR.
///
/// Example
/// -------
///
/// - Create geometry description for 2 rings of 48 detectors using 1 million
///   emissions from each pixel:
///
///        2d_barrel_geometry -s square -w 0.007 -h 0.017
///                           -r 0.360 -d 48 --radius2 0.400
///                           -e 1000000 -o data/201412_rings/gpu_2rings
///
/// Authors
/// -------
/// - Piotr Bialas <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/geometry_cmd.txt
///
/// \sa \ref cmd_2d_barrel_phantom, \ref cmd_2d_barrel_lm_reconstruction

#include <iostream>
#include <fstream>
#include <deque>
#include <random>

#include "2d/barrel/boost_geometry_utils.h"
#include "2d/barrel/options.h"

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/progress.h"
#include "util/backtrace.h"

#include "scanner_builder.h"
#include "ring_scanner.h"
#include "generic_scanner.h"

#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/barrel/geometry.h"

using F = float;
using S = short;
using RNG = std::mt19937;

using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Point = PET2D::Point<F>;

using Geometry = PET2D::Barrel::Geometry<F, S>;
using LORInfo = PET2D::Barrel::LORInfo<F, S>;
using PixelInfo = PET2D::Barrel::Geometry<F, S>::PixelInfo;
using PixelInfoContainer = PET2D::Barrel::LORInfo<F, S>::PixelInfoList;
using LOR = PET2D::Barrel::LOR<S>;

using BoostGeometryUtils = PET2D::Barrel::BoostGeometryUtils<F, S>;

using Polygon = typename BoostGeometryUtils::Polygon;
using Point2D = BoostGeometryUtils::Point2D;

int main(int argc, char* argv[]) {
  try {
    cmdline::parser cl;
    PET2D::Barrel::add_matrix_options(cl);
    cl.try_parse(argc, argv);

    // FIXME: this only works for big barrel
    PET2D::Barrel::set_big_barrel_options(cl);
    auto scanner =
        PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
            PET2D_BARREL_SCANNER_CL(cl, F));

    auto verbose = cl.exist("verbose");
    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();

    std::vector<Polygon> detectors;
    std::vector<Point> detectors_centers;

    F pixel_size = F(0.005);
    if (cl.exist("s-pixel"))
      pixel_size = cl.get<double>("s-pixel");
    F fov_radius = cl.get<double>("fov-radius");

    if (verbose) {
      std::cout << "fov " << fov_radius << " size " << pixel_size << std::endl;
    }

    S n_columns, n_rows;
    if (!cl.exist("n-pixels")) {
      n_columns = 2 * S(std::ceil(fov_radius / pixel_size));
    } else {
      n_columns = cl.get<int>("n-pixels");
    }
    n_rows = n_columns;

    if (verbose) {
      std::cout << "cols " << n_columns << " rows " << n_rows << std::endl;
    }

    PET2D::PixelGrid<F, S> grid(
        n_columns,
        n_rows,
        pixel_size,
        Point(-pixel_size * n_columns / 2, -pixel_size * n_rows / 2));

    for (int i = 0; i < (int)scanner.size(); i++) {
      auto detector = scanner[i];
      Polygon detector_poly = BoostGeometryUtils::make_detector(detector);

      detectors.push_back(detector_poly);
      detectors_centers.push_back(detector.center());
    }

    std::ofstream svg(output_base_name + "_map.svg");
    boost::geometry::svg_mapper<Point2D> mapper(svg, 1200, 1200);

    for (const auto& detector : detectors) {
#if DEBUG
      std::cout << boost::geometry::wkt(detector) << std::endl;
#endif
      mapper.add(detector);
    }

    auto fov_circle = BoostGeometryUtils::make_circle(
        Point(0, 0), cl.get<double>("fov-radius"), 128);
    mapper.add(fov_circle);

    for (const auto& detector : detectors) {
      mapper.map(detector, "fill:rgb(0,0,255);");
    }
    mapper.map(fov_circle, "fill:none;stroke:red;");

    S n_detectors = scanner.size();

    std::vector<LOR> lor_map;
    auto lor_count = LOR::end_for_detectors(n_detectors).index();
    lor_map.resize(lor_count);
    for (LOR lor(0, 0); lor < LOR::end_for_detectors(n_detectors); ++lor) {
      lor_map[lor.index()] = lor;
    }

    Geometry geometry(n_detectors, grid);

    util::progress progress(verbose, lor_count);

#if _OPENMP && !_MSC_VER
// We need to try catch inside OpenMP thread, otherwise we will not see the
// error thrown.
#define TRY try {
#define CATCH                     \
  }                               \
  catch (std::string & ex) {      \
    std::cerr << ex << std::endl; \
    throw(ex);                    \
  }
#else
#define TRY
#define CATCH
#endif

#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int lor_index = 0; lor_index < lor_count; ++lor_index) {
      auto lor = lor_map[lor_index];
      progress(lor_index);
      TRY;
      boost::geometry::model::multi_polygon<Polygon> pair;

      boost::geometry::union_(detectors[lor.first],   // combine
                              detectors[lor.second],  // these
                              pair);
      Polygon lor_hull;
      boost::geometry::convex_hull(pair, lor_hull);

      PET2D::LineSegment<F> segment(detectors_centers[lor.second],
                                    detectors_centers[lor.first]);

      // TODO: Calculate width of the LOR.
      auto width1 = F(0);
      auto width2 = F(0);
      Detector detector1 = scanner[lor.first];
      Detector detector2 = scanner[lor.second];
      for (int i = 0; i < (int)detector1.size(); ++i) {
        auto p1 = detector1[i];
        auto dist1 = std::abs(segment.distance_from(p1));
        if (dist1 > width1)
          width1 = dist1;

        auto p2 = detector2[i];
        auto dist2 = std::abs(segment.distance_from(p2));
        if (dist2 > width2)
          width2 = dist2;
      }

      LORInfo lor_info(lor, segment, width1 + width2);

      if (boost::geometry::intersects(lor_hull, fov_circle)) {
        for (int ix = 0; ix < grid.n_columns; ++ix)
          for (int iy = 0; iy < grid.n_rows; ++iy) {

            Point center = grid.center_at(ix, iy);
            Polygon pixel = BoostGeometryUtils::make_pixel(grid, ix, iy);
            if (boost::geometry::intersects(pixel, fov_circle)) {
              boost::geometry::model::multi_polygon<Polygon> inter;
              boost::geometry::intersection(lor_hull, pixel, inter);
              auto area = boost::geometry::area(inter);

              if (area > 0) {
                auto pixel_area = boost::geometry::area(pixel);
                auto fill = area / pixel_area;
                auto t = segment.projection_scaled(center);
                auto distance = segment.distance_from(center);
                PixelInfo pixel_info;
                pixel_info.pixel = PET2D::Pixel<S>(ix, iy);
                pixel_info.t = t;
                pixel_info.distance = distance;
                pixel_info.fill = fill;
                lor_info.push_back(pixel_info);
              }
            }
          }

        lor_info.sort();
      }
      geometry[lor] = lor_info;

      CATCH;
      progress(lor_index, true);
    }

    util::obstream out_geometry(output);
    out_geometry << geometry;

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  }
  return 1;
}

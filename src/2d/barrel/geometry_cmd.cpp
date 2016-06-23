/// \page cmd_2d_barrel_geometry 2d_barrel_geometry
/// \brief 2D Barrel PET geometry description construction tool
///
/// Creates geometry description binary file for LM reconstruction
/// \ref cmd_2d_barrel_lm_reconstruction.
///
/// \warning This command is experimental and should be **NOT** used for regular
/// 2D bin-mode reconstruction.
///
/// This is alternative for \ref cmd_2d_barrel_matrix. It does not use
/// Monte-Carlo, but calculates every LOR geometry and pixels that lie inside
/// this LOR.
///
/// Example
/// -------
///
/// 1. Make a \c playground directory and step into it
///
///        mkdir playground
///        cd playground
///
/// 2. Create geometry description for 2 rings of 48 detectors using 1 million
///    emissions from each pixel:
///
///        ../2d_barrel_geometry \
///          -s square -w 0.007 -h 0.017 \
///          -r 0.360,0.400 -d 48 \
///          -e 1m \
///          -o g_2rings
///
/// Authors
/// -------
/// - Piotr Bialas <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verboutput 2d_barrel_geometry
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
#include "3d/hybrid/options.h"

#include "common/types.h"

using RNG = std::mt19937;

using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;

using Geometry = PET2D::Barrel::Geometry<F, S>;
using LORGeometry = PET2D::Barrel::LORGeometry<F, S>;
using PixelInfo = PET2D::Barrel::Geometry<F, S>::PixelInfo;
using PixelInfoContainer = PET2D::Barrel::LORGeometry<F, S>::PixelInfoList;
using LOR = PET2D::Barrel::LOR<S>;

using BoostGeometryUtils = PET2D::Barrel::BoostGeometryUtils<F, S>;

using Polygon = typename BoostGeometryUtils::Polygon;
using Point2D = BoostGeometryUtils::Point2D;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  cl.add<double>("length", 0, " barrel length - not used !!", false, 0.0);
  cl.add<double>(
      "z-position", 'z', "z plane position - not used !!", false, 0.0);
  cl.parse_check(argc, argv);
  PET2D::Barrel::calculate_scanner_options(cl, argc);

  auto scanner = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, F));

  auto verbose = cl.count("verbose");
  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();

  std::vector<Polygon> detectors;
  std::vector<Point> detectors_centers;

  auto pixel_size = cl.get<double>("s-pixel");
  auto fov_radius = cl.get<double>("fov-radius");

  if (verbose) {
    std::cout << " fov: " << fov_radius << std::endl
              << "size: " << pixel_size << std::endl;
  }

  S n_columns, n_rows;
  if (!cl.exist("n-pixels")) {
    n_columns = 2 * S(std::ceil(fov_radius / pixel_size));
  } else {
    n_columns = cl.get<int>("n-pixels");
  }
  n_rows = n_columns;

  if (verbose) {
    std::cout << "cols: " << n_columns << std::endl
              << "rows: " << n_rows << std::endl;
  }

  PET2D::PixelGrid<F, S> grid(n_columns, n_rows, pixel_size);

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

    LORGeometry lor_geometry(lor, segment, width1 + width2);

    if (boost::geometry::intersects(lor_hull, fov_circle)) {
      for (int ix = 0; ix < grid.n_columns; ++ix)
        for (int iy = 0; iy < grid.n_rows; ++iy) {
          Pixel pixel_coord(ix, iy);
          Point center = grid.center_at(pixel_coord);
          Polygon pixel = BoostGeometryUtils::make_pixel(grid, pixel_coord);
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
              pixel_info.weight = F(0);
              lor_geometry.push_back(pixel_info);
            }
          }
        }

      lor_geometry.sort();
    }
    geometry[lor] = lor_geometry;

    CATCH;
    progress(lor_index, true);
  }

  util::obstream out_geometry(output);
  out_geometry << geometry;

  CMDLINE_CATCH
}

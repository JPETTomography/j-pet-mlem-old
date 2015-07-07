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
#include "scanner_builder.h"
#include "ring_scanner.h"
#include "generic_scanner.h"
#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/barrel/lors_pixels_info.h"

using F = float;
using S = short;
using RNG = std::mt19937;

using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Point = PET2D::Point<F>;

using PixelInfo = PET2D::Barrel::LORsPixelsInfo<F, S>::PixelInfo;
using PixelInfoContainer =
    PET2D::Barrel::LORsPixelsInfo<F, S>::PixelInfoContainer;
using LOR = PET2D::Barrel::LOR<S>;

using BoostGeometryUtils = PET2D::Barrel::BoostGeometryUtils<F, S>;

using Polygon = typename BoostGeometryUtils::Polygon;
using Point2D = BoostGeometryUtils::Point2D;

int main(int argc, char* argv[]) {
  try {
    cmdline::parser cl;
    PET2D::Barrel::add_matrix_options(cl);
    cl.try_parse(argc, argv);

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

    std::cout << "fov " << fov_radius << " size " << pixel_size << "\n";
    S n_columns, n_rows;
    if (!cl.exist("n-pixels")) {
      n_columns = 2 * S(std::ceil(fov_radius / pixel_size));
    } else {
      n_columns = cl.get<int>("n-pixels");
    }
    n_rows = n_columns;
    std::cout << "cols " << n_columns << " rows " << n_rows << "\n";

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
    PET2D::Barrel::LORsPixelsInfo<F, S> lor_info(n_detectors, grid);

    std::ofstream lor_info_stream(output, std::ios::binary);
    lor_info_stream.write((const char*)&n_detectors, sizeof(S));
    grid.write(lor_info_stream);

    std::vector<LOR> lor_map;
    auto lor_count = LOR::end_for_detectors(n_detectors).index();
    lor_map.resize(lor_count);
    for (LOR lor(0, 0); lor < LOR::end_for_detectors(n_detectors); ++lor) {
      lor_map[lor.index()] = lor;
    }

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

      boost::geometry::union_(
          detectors[lor.first], detectors[lor.second], pair);
      Polygon lor_polygon;
      boost::geometry::convex_hull(pair, lor_polygon);

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
      F width = width1 + width2;

      if (boost::geometry::intersects(lor_polygon, fov_circle)) {
        for (int ix = 0; ix < grid.n_columns; ++ix)
          for (int iy = 0; iy < grid.n_rows; ++iy) {

            Point center = grid.center_at(ix, iy);
            Polygon pixel = BoostGeometryUtils::make_pixel(grid, ix, iy);
            if (boost::geometry::intersects(pixel, fov_circle)) {
              boost::geometry::model::multi_polygon<Polygon> inter;
              boost::geometry::intersection(lor_polygon, pixel, inter);
              auto area = boost::geometry::area(inter);

              if (area > 0) {
                auto pixel_area = boost::geometry::area(pixel);
                auto fill = area / pixel_area;
                auto t = segment.projection_scaled(center);
                auto distance = segment.distance_from(center);
                PixelInfo info;
                info.pixel = PET2D::Pixel<S>(ix, iy);
                info.t = t;
                info.distance = distance;
                info.fill = fill;
                lor_info.push_back_pixel_info(lor, info);
              }
            }
          }

        lor_info.sort();
      }

      int n_pixels = lor_info[lor].pixels.size();

#if _OPENMP
#pragma omp critical
#endif
      {
        lor_info_stream.write((const char*)&lor.first, sizeof(lor.first));
        lor_info_stream.write((const char*)&lor.second, sizeof(lor.second));
        lor_info_stream.write((const char*)&detectors_centers[lor.first],
                              2 * sizeof(F));
        lor_info_stream.write((const char*)&detectors_centers[lor.second],
                              2 * sizeof(F));
        lor_info_stream.write((const char*)&width, sizeof(F));
        lor_info_stream.write((const char*)&n_pixels, sizeof(int));
        if (n_pixels > 0)
          lor_info_stream.write((const char*)&lor_info[lor].pixels[0],
                                n_pixels * sizeof(PixelInfo));
      }
      CATCH;
      progress(lor_index, true);
    }

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
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

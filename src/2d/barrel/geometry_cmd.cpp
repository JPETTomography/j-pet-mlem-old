
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

    for (int i = 0; i < scanner.size(); i++) {
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

    // Loop over the LORs
    int count = 0;
    util::progress progress(verbose, n_detectors * (n_detectors + 1) / 2);
    for (S d1 = 0; d1 < n_detectors; ++d1) {
      for (S d2 = 0; d2 < d1; ++d2, ++count) {
        progress(count);

        boost::geometry::model::multi_polygon<Polygon> pair;

        boost::geometry::union_(detectors[d1], detectors[d2], pair);
        Polygon lor;
        boost::geometry::convex_hull(pair, lor);

#if DEBUG
        std::cout << "l: " << count << "  " << d1 << " " << d2 << "\n";
#endif
        lor_info_stream.write((const char*)&d1, sizeof(d1));
        lor_info_stream.write((const char*)&d2, sizeof(d2));

        PET2D::LineSegment<F> segment(detectors_centers[d2],
                                      detectors_centers[d1]);

        // TODO: Calculate width of the LOR.
        auto width1 = F(0);
        auto width2 = F(0);
        Detector detector1 = scanner[d1];
        Detector detector2 = scanner[d2];
        for (int i = 0; i < detector1.size(); ++i) {
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
        lor_info_stream.write((const char*)&detectors_centers[d1],
                              2 * sizeof(F));
        lor_info_stream.write((const char*)&detectors_centers[d2],
                              2 * sizeof(F));
#if DEBUG
        std::cout << detectors_centers[d1].x << " " << detectors_centers[d1].y
                  << " "  //
                  << detectors_centers[d2].x << " " << detectors_centers[d2].y
                  << "\n";
#endif
        lor_info_stream.write((const char*)&width, sizeof(F));
        if (boost::geometry::intersects(lor, fov_circle)) {
          for (int ix = 0; ix < grid.n_columns; ++ix)
            for (int iy = 0; iy < grid.n_rows; ++iy) {

              Point center = grid.center_at(ix, iy);
              Polygon pixel = BoostGeometryUtils::make_pixel(grid, ix, iy);
              if (boost::geometry::intersects(pixel, fov_circle)) {
                boost::geometry::model::multi_polygon<Polygon> inter;
                boost::geometry::intersection(lor, pixel, inter);
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
                  lor_info.push_back_pixel_info(LOR(d1, d2), info);
                }
              }
            }

          lor_info.sort();
        }

        int n_pixels = lor_info[LOR(d1, d2)].pixels.size();

        lor_info_stream.write((const char*)&n_pixels, sizeof(int));
        if (n_pixels > 0)
          lor_info_stream.write((const char*)&lor_info[LOR(d1, d2)].pixels[0],
                                n_pixels * sizeof(PixelInfo));

        progress(count, true);
      }
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

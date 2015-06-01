
#include <iostream>
#include <deque>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/geometries.hpp>
//#include <boost/geometry/strategies/agnostic/hull_graham_andrew.hpp>

#include <boost/foreach.hpp>

#include "2d/barrel/options.h"

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "scanner_builder.h"
#include "ring_scanner.h"
#include "generic_scanner.h"
#include "2d/geometry/line_segment.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/barrel/lor_info.h"

using FType = float;
using SType = int;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 192, SType>;
using Point = PET2D::Point<FType>;
using point_2d = boost::geometry::model::d2::point_xy<FType>;

using Polygon = boost::geometry::model::polygon<point_2d>;

Polygon makeCircle(const Point& center, FType radius, int n = 64) {
  Polygon circle;
  FType da = 2 * M_PI / n;
  FType angle = 0.0;
  for (int i = 0; i < n; ++i) {
    FType x, y;
    x = radius * std::cos(angle);
    y = radius * std::sin(angle);
    boost::geometry::append(circle, boost::geometry::make<point_2d>(x, y));
    angle += da;
  }
  boost::geometry::append(circle,
                          boost::geometry::make<point_2d>(radius, FType(0.0)));
  return circle;
}

Polygon makePixel(const PET2D::PixelGrid<FType, SType>& grid, int ix, int iy) {
  Polygon pixel;
  auto size = grid.pixel_size;
  Point ll = grid.lower_left_at(ix, iy);
  auto x = ll.x;
  auto y = ll.y;
  // std::cout<<x<<" "<<y<<"\n";
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y + size));
  boost::geometry::append(pixel,
                          boost::geometry::make<point_2d>(x + size, y + size));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x + size, y));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
  return pixel;
}

struct PixelInfo {
  PET2D::Pixel<SType> pixel;
  FType t;
  FType distance;
  FType fill;
};

int main(int argc, char* argv[]) {
  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  try {
    cl.try_parse(argc, argv);
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

  PET2D::Barrel::set_big_barrel_options(cl);
  auto scanner = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, FType));

  std::vector<Polygon> detectors;
  std::vector<Point> detectors_centers;

  FType pixel_size = 0.005;
  if (cl.exist("s-pixel"))
    pixel_size = cl.get<double>("s-pixel");
  FType fov_radius = cl.get<double>("fov-radius");
  std::cout << "fov " << fov_radius << " size " << pixel_size << "\n";
  SType n_columns = 2 * SType(std::ceil(fov_radius / pixel_size));
  SType n_rows = n_columns;
  std::cout << "cols " << n_columns << " rows " << n_rows << "\n";
  PET2D::PixelGrid<FType, SType> grid(
      n_columns,
      n_rows,
      pixel_size,
      Point(-pixel_size * n_columns / 2, -pixel_size * n_rows / 2));

  for (int i = 0; i < scanner.size(); i++) {
    auto detector = scanner[i];
    Polygon detector_poly;
    for (int j = 0; j < detector.size(); j++) {
      Point p = detector[j];
      boost::geometry::append(detector_poly,
                              boost::geometry::make<point_2d>(p.x, p.y));
    }
    Point p = detector[0];
    boost::geometry::append(detector_poly,
                            boost::geometry::make<point_2d>(p.x, p.y));
    detectors.push_back(detector_poly);
    detectors_centers.push_back(detector.center());
  }

  std::ofstream svg("my_map.svg");
  boost::geometry::svg_mapper<point_2d> mapper(svg, 1200, 1200);

  for (auto p = detectors.begin(); p != detectors.end(); ++p) {
    // std::cout << boost::geometry::wkt(*p) << std::endl;
    mapper.add(*p);
  }

  auto fov_circle = makeCircle(Point(0, 0), cl.get<double>("fov-radius"), 128);
  mapper.add(fov_circle);

  for (auto p = detectors.begin(); p != detectors.end(); ++p) {
    mapper.map(*p, "fill:rgb(0,0,255);");
  }
  mapper.map(fov_circle, "fill:none;stroke:red;");
  //mapper.~svg_mapper();

  int n_detectors = scanner.size();
  int i = 0;
  for (int d1 = 0; d1 < n_detectors; ++d1) {
    for (int d2 = 0; d2 < d1; ++d2) {
      boost::geometry::model::multi_polygon<Polygon> pair;

      boost::geometry::union_(detectors[d1], detectors[d2], pair);
      Polygon lor;
      boost::geometry::convex_hull(pair, lor);

      if (boost::geometry::intersects(lor, fov_circle)) {
        std::cout << "l : " << i << "  " << d1 << " " << d2 << "\n";
        i++;
        std::vector<PixelInfo> pixel_info;
        PET2D::LineSegment<FType> segment(detectors_centers[d2],
                                          detectors_centers[d1]);

        for (int ix = 0; ix < grid.n_columns; ++ix)
          for (int iy = 0; iy < grid.n_rows; ++iy) {

            Point center = grid.center_at(ix, iy);
            Polygon pixel = makePixel(grid, ix, iy);
            if (boost::geometry::intersects(pixel, fov_circle)) {
              boost::geometry::model::multi_polygon<Polygon> inter;
              boost::geometry::intersection(lor, pixel, inter);
              auto area = boost::geometry::area(inter);

              if (area > 0) {
                auto pixel_area = boost::geometry::area(pixel);
                auto fill = area / pixel_area;
                auto t = segment.projection_on(center);
                auto distance = segment.distance_from(center);
                PixelInfo info;
                info.pixel = PET2D::Pixel<SType>(ix, iy);
                info.t = t;
                info.distance = distance;
                info.fill = fill;
                pixel_info.push_back(info);
              }
            }
          }

        std::sort(
            pixel_info.begin(),
            pixel_info.end(),
            [](const PixelInfo& a, const PixelInfo& b) { return a.t < b.t; });

        for (PixelInfo info : pixel_info) {
          std::cout << "p : " << info.pixel.x << " " << info.pixel.y << " "
                    << info.t << info.fill << " "
                    << " " << info.distance << "\n";
        }
      }
    }
  }
}

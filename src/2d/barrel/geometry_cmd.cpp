
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

using FType = float;
using SType = short;
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

Polygon makePixel(FType size, int ix, int iy) {
  Polygon pixel;
  FType x = ix * size;
  FType y = iy * size;
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y + size));
  boost::geometry::append(pixel,
                          boost::geometry::make<point_2d>(x + size, y + size));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x + size, y));
  boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
  return pixel;
}

int main(int argc, char* argv[]) {
  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  cl.try_parse(argc, argv);

  PET2D::Barrel::set_big_barrel_options(cl);
  auto scanner = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, FType));

  std::vector<Polygon> detectors;

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

  int n_detectors = 192;
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
        for (int ix = -50; ix < 50; ++ix)
          for (int iy = -50; iy < 50; ++iy) {
            Polygon pixel = makePixel(0.005, ix, iy);
            if (boost::geometry::intersects(pixel, fov_circle)) {
              boost::geometry::model::multi_polygon<Polygon> inter;
              boost::geometry::intersection(lor, pixel, inter);
              auto area = boost::geometry::area(inter);
              auto pixel_area = boost::geometry::area(pixel);
              auto fill = area / pixel_area;
              if (area > 0) {
                std::cout << "p : " << ix << " " << iy << " " << fill << "\n";
                //                mapper.add(pixel);
                //                mapper.add(lor);

                //                auto fill = (int)floor(255*area/pixel_area);
                //                char fill_style[64];
                //                sprintf(fill_style,"fill:rgb(%d,255, %d);",
                //                255-fill, 255-fill);
                //                mapper.map(pixel, fill_style);
              }
            }
          }
        //mapper.map(lor, "fill:none;stroke:rgb(0,0,255);");
        // goto end;
      }
    }
  }
end:
  ;

  mapper.map(fov_circle, "fill:none;stroke:red;");
}

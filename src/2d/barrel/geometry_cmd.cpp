
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

int main(int argc, char* argv[]) {
  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  cl.try_parse(argc, argv);

  PET2D::Barrel::set_big_barrel_options(cl);
  auto scanner = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, FType));

  typedef boost::geometry::model::polygon<point_2d> Polygon;

  std::vector<Polygon> detectors;

  for (int i = 0; i < scanner.size(); i++) {
    auto detector = scanner[i];
    Polygon detector_poly;
    for (int j = 0; j < detector.size(); j++) {
      Point p = detector[j];
      boost::geometry::append(detector_poly,
                              boost::geometry::make<point_2d>(p.x, p.y));
    }
    detectors.push_back(detector_poly);
  }

  std::ofstream svg("my_map.svg");
  boost::geometry::svg_mapper<point_2d> mapper(svg, 1024, 1024);

  for (auto p = detectors.begin(); p != detectors.end(); ++p) {
    std::cout << boost::geometry::wkt(*p) << std::endl;
    mapper.add(*p);
  }

  boost::geometry::model::multi_polygon<Polygon> pair;

  boost::geometry::union_(detectors[0], detectors[20],pair);
  Polygon lor;
  // boost::geometry::convex_hull(pair, lor);
  mapper.add(pair);
  mapper.map(pair, "fill::rgb(1,0,0)");
  std::cout << boost::geometry::wkt(pair) << std::endl;
  for (auto p = detectors.begin(); p != detectors.end(); ++p) {
     mapper.map(*p, "fill::rgb(0,0,0)");
  }


}

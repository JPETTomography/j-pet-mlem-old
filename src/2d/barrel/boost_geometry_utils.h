#pragma once

#ifdef HAVE_Boost

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include "2d/geometry/pixel_grid.h"
namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class BoostGeometryUtils {
 public:
  using point_2d = boost::geometry::model::d2::point_xy<FType>;
  using Polygon = boost::geometry::model::polygon<point_2d>;
  using Point = PET2D::Point<FType>;

  static Polygon makePixel(const PET2D::PixelGrid<FType, SType>& grid,
                           int ix,
                           int iy) {
    Polygon pixel;
    auto size = grid.pixel_size;
    Point ll = grid.lower_left_at(ix, iy);
    auto x = ll.x;
    auto y = ll.y;
    // std::cout<<x<<" "<<y<<"\n";
    boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
    boost::geometry::append(pixel,
                            boost::geometry::make<point_2d>(x, y + size));
    boost::geometry::append(
        pixel, boost::geometry::make<point_2d>(x + size, y + size));
    boost::geometry::append(pixel,
                            boost::geometry::make<point_2d>(x + size, y));
    boost::geometry::append(pixel, boost::geometry::make<point_2d>(x, y));
    return pixel;
  }

  static Polygon makeCircle(const Point& center, FType radius, int n = 64) {
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
    boost::geometry::append(
        circle, boost::geometry::make<point_2d>(radius, FType(0.0)));
    return circle;
  }

  template <typename Detector>
  static Polygon makeDetector(const Detector& detector) {
    Polygon detector_poly;
    for (int j = 0; j < detector.size(); j++) {
      Point p = detector[j];
      boost::geometry::append(detector_poly,
                              boost::geometry::make<point_2d>(p.x, p.y));
    }
    Point p = detector[0];
    boost::geometry::append(detector_poly,
                            boost::geometry::make<point_2d>(p.x, p.y));
    return detector_poly;
  }

  static Polygon makeLor(const Polygon& d1, const Polygon& d2) {
    boost::geometry::model::multi_polygon<Polygon> pair;

    boost::geometry::union_(d1, d2, pair);
    Polygon lor;
    boost::geometry::convex_hull(pair, lor);
    return lor;
  }

  template <typename Detector>
  static Polygon makeLor(const Detector& d1, const Detector& d2) {
    return makeLor(makeDetector(d1), makeDetector(d2));
  }
};
}
}
#endif

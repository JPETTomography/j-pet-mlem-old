#pragma once

#ifdef HAVE_BOOST

#ifdef __GNUC__
#if __GNUC__ > 3 && __GNUC_MINOR__ > 6
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-local-typedef"
#endif
#endif
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wunused-local-typedef"
#endif

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/geometries.hpp>

#include <boost/foreach.hpp>

#ifdef __GNUC__
#if __GNUC__ > 3 && __GNUC_MINOR__ > 6
#pragma GCC diagnostic pop
#endif
#endif
#ifdef __clang__
#pragma clang diagnostic pop
#endif

#include "2d/geometry/pixel_grid.h"

namespace PET2D {
namespace Barrel {

template <typename FType, typename SType> class BoostGeometryUtils {
 public:
  using F = FType;
  using S = SType;
  using point_2d = boost::geometry::model::d2::point_xy<F>;
  using Polygon = boost::geometry::model::polygon<point_2d>;
  using Point = PET2D::Point<F>;

  static Polygon make_pixel(const PET2D::PixelGrid<F, S>& grid, S ix, S iy) {
    Polygon pixel;
    auto size = grid.pixel_size;
    Point ll = grid.lower_left_at(ix, iy);
    auto x = ll.x;
    auto y = ll.y;

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

  static Polygon make_circle(const Point& center, F radius, int n = 64) {
    Polygon circle;
    F da = 2 * M_PI / n;
    F angle = 0.0;
    for (int i = 0; i < n; ++i) {
      F x, y;
      x = center.x + radius * std::cos(angle);
      y = center.y + radius * std::sin(angle);
      boost::geometry::append(circle, boost::geometry::make<point_2d>(x, y));
      angle += da;
    }
    boost::geometry::append(circle,
                            boost::geometry::make<point_2d>(radius, F(0)));
    return circle;
  }

  template <typename Detector>
  static Polygon make_detector(const Detector& detector) {
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

  static Polygon make_lor(const Polygon& d1, const Polygon& d2) {
    boost::geometry::model::multi_polygon<Polygon> pair;

    boost::geometry::union_(d1, d2, pair);
    Polygon lor;
    boost::geometry::convex_hull(pair, lor);
    return lor;
  }

  template <typename Detector>
  static Polygon make_lor(const Detector& d1, const Detector& d2) {
    return make_lor(make_detector(d1), make_detector(d2));
  }
};
}
}
#endif

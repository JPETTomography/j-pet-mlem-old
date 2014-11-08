#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/circle.h"
#include "util/array.h"
#include "util/cuda/compat.h"
#if !__CUDACC__
#include "util/svg_ostream.h"
#endif

namespace PET2D {
namespace Barrel {

/// Circular shape detector
template <typename FType = double> class CircleDetector : public Circle<FType> {
 public:
  using F = FType;
  using Base = Circle<F>;
  using Angle = F;
  using Point = PET2D::Point<F>;
  using Intersections = util::array<2, Point>;
  using Event = typename Base::Event;

  CircleDetector() = delete;

  CircleDetector(F radius) : Circle<F>(radius), center(0, 0) {}

  // this is for compatibility with square detector
  CircleDetector(F w, F h, F d) : Base(w / 2), center(0, 0) {
    (void)d;  // unused
    if (w != h)
      throw("circle detector height and width must be equal");
  }

  static F default_height_for_width(const F w) { return w; }

  CircleDetector(F radius, const Point& center)
      : Base(radius), center(center) {}

  CircleDetector& rotate(Angle phi) {
    center.rotate(phi);
    return *this;
  }

  CircleDetector& operator+=(Point t) {
    center += t;
    return *this;
  }

  F max_distance() { return center.length() + this->radius; }

  Point center;

  /// \returns itself
  const CircleDetector& circumscribe_circle() const { return *this; }

  _ Intersections intersections(typename Base::Event e) {
    auto intersections = this->secant(e - center);
    for (auto& p : intersections) {
      p += center;
    }
    return intersections;
  }

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          CircleDetector& cd) {
    svg << "<circle cx=\"" << cd.center.x << "\" cy=\"" << cd.center.y
        << "\" r=\"" << cd.radius << "\"/>" << std::endl;
    return svg;
  }

  friend std::ostream& operator<<(std::ostream& out, CircleDetector& cd) {
    out << "circle (" << cd.center.x << ", " << cd.center.y << ") radius "
        << cd.radius << std::endl;
    out << std::flush;
    return out;
  }
#endif
};
}  // Barrel
}  // PET2D

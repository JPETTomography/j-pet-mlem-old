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

/// Represents single detector with round shape:
/// \image html shape_circle.pdf.png
template <typename FType = double> class CircleDetector : public Circle<FType> {
 public:
  using F = FType;
  using Base = Circle<F>;
  using Angle = F;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Intersections = util::array<2, Point>;
  using Event = typename Base::Event;

  CircleDetector() = delete;

  CircleDetector(F radius) : Base(radius), center(0, 0) {}

  // this is for compatibility with square detector
  CircleDetector(F w, F h, F d) : Base(w / 2), center(0, 0) {
    (void)d;  // unused
    if (w != h)
      throw("circle detector height and width must be equal");
  }

  static F default_height_for_width(const F w) { return w; }

  CircleDetector(F radius, const Point& center_a)
      : Base(radius), center(center_a) {}

  CircleDetector& rotate(Angle phi) {
    center.rotate(phi);
    return *this;
  }

  CircleDetector& operator+=(Vector t) {
    center += (t);
    return *this;
  }

  CircleDetector operator+(Vector t) const {
    CircleDetector result(*this);
    result.center += t;
    return reinterpret_cast<CircleDetector&&>(result);
  }

  F max_distance() { return center.distance_from_origin() + this->radius; }

  /// \returns itself
  const CircleDetector& circumscribe_circle() const { return *this; }

  _ bool intersects(typename Base::Event e) const {
    return Base::intersects(e - center.as_vector());
  }

  _ Intersections intersections(typename Base::Event e) const {
    Vector v = center.as_vector();
    auto intersections = this->secant(e - v);
    for (auto& p : intersections) {
      p += v;
    }
    return intersections;
  }

 public:
  Point center;

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          CircleDetector& cd) {
    svg << "<circle cx=\"" << cd.center.x << "\" cy=\"" << cd.center.y
        << "\" r=\"" << cd.radius << "\"/>" << std::endl;
    return svg;
  }

  friend std::ostream& operator<<(std::ostream& out, CircleDetector& cd) {
    out << "circle (" << cd.x << ", " << cd.y << ") radius " << cd.radius
        << std::endl;
    out << std::flush;
    return out;
  }
#endif
};
}  // Barrel
}  // PET2D

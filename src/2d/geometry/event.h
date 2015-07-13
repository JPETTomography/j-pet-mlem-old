#pragma once

#if !__CUDACC__
#include <iostream>
#endif

#include "2d/geometry/point.h"
#include "2d/geometry/vector.h"

namespace PET2D {

/// Generic emitted 2D event

/// Events of this type  are emmited by every 2D phantom.
template <typename FType> struct Event {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;

  _ Event(const Point& center, const Vector& direction)
      : center(center), direction(direction) {}

  _ Event(F x, F y, F dx, F dy) : Event(Point(x, y), Vector(dx, dy)) {}

  _ Event(const Point& center, F angle)
      : Event(center, Vector(compat::sin(angle), compat::cos(angle))) {}

  const Point center;
  const Vector direction;

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Event& event) {
    out << event.center << " " << event.direction;
    return out;
  }
#endif
};

}  // PET2D

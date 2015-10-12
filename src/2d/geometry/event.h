#pragma once

#if !__CUDACC__
#include <iostream>
#endif

#include "2d/geometry/point.h"
#include "2d/geometry/vector.h"

namespace PET2D {

/// Generic 2D emission event
////
/// Emission event consist of origin (emission) point \f$ origin = (x, y) \f$
/// and direction (vector) \f$ direction = (dx, dy) \f$.
///
/// Direction vector can be constructed from \f$ \phi \f$ angle using:
/// \f[
///   direction = (sin(\phi), cos(\phi))
/// \f]
template <typename FType> struct Event {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;

  /// Emission event at \f$ origin = (x, y) \f$ point with
  /// \f$ direction = (dx, dy) \f$ direction.
  _ Event(const Point& origin, const Vector& direction)
      : origin(origin), direction(direction) {}

  /// Emission event at \f$ (x, y) \f$ point with \f$ (dx, dy) \f$ direction.
  _ Event(F x, F y, F dx, F dy) : Event(Point(x, y), Vector(dx, dy)) {}

  /// Emission event at \f$ origin = (x, y) \f$ point and \f$ \phi \f$ angle.
  ////
  /// \note Angle is counted from \f$ (0, 1) \f$ and grows clockwise.
  /// \todo FIXME: This is somehow different than PET2D::Barrel::Event where
  /// angle counted from \f$ (1, 0) \f$ and grows counter-clockwise.
  _ Event(const Point& origin, F phi)
      : Event(origin, Vector(compat::sin(phi), compat::cos(phi))) {}

  const Point origin;
  const Vector direction;

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Event& event) {
    out << event.origin << " " << event.direction;
    return out;
  }
#endif
};

}  // PET2D

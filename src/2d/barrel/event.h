#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/vector.h"
#include "2d/geometry/event.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// 2D barrel PET emission event

/// This is extension to the PET2D::Event generic 2D event holding normal to the
/// direction and distance to origin general line equation
/// \f$ a x + b y + c = 1 \f$ coefficients, and also some pre-computed forms.
///
/// \note We do not normalize \f$ direction = (dx, dy)\f$ vector, since it
/// comes usually from the angle, so \f$ normal = (a, b) = (-dy, dx) =
/// sin(\phi), b = -cos(\phi) \f$ then \f$ a^2 + b^2 = 1 \f$.
template <typename FType> struct Event : public PET2D::Event<FType> {
  using F = FType;

  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Base = PET2D::Event<FType>;

  /// Event requires usage of a concrete constructor.
  Event() = delete;

  _ Event(const Point& origin, const Vector& direction)
      : Base(origin, direction),
        normal(direction.cw_perpendicular()),
        // line equation c coefficient: a x + b y == c
        distance_from_origin(origin.as_vector().dot(normal)),
        // helper variables
        b2c(normal.y * normal.y * distance_from_origin),
        ac(normal.x * distance_from_origin),
        c2(distance_from_origin * distance_from_origin),
        inv_b(1 / normal.y) {}

  _ Event(const Base& event) : Event(event.origin, event.direction) {}

  /// Emission event at \f$ origin = (x, y) \f$ point and \f$ \phi \f$ angle.

  /// \note Angle is counted from \f$ (1, 0) \f$ and grows counter-clockwise.
  /// \todo FIXME: This is somehow different than PET2D::Event where
  /// angle counted from \f$ (0, 1) \f$ and grows clockwise.
  _ Event(const Point origin, F phi)
      : Event(origin, Vector(std::cos(phi), std::sin(phi))) {}

  /// Emission event at \f$ (x, y) \f$ point and \f$ \phi \f$ angle.

  /// \note Angle is counted from \f$ (1, 0) \f$ and grows counter-clockwise.
  _ Event(F x, F y, F phi) : Event(Point(x, y), phi) {}

  /// Emission event at \f$ (x, y) \f$ point with \f$ (dx, dy) \f$ direction.
  _ Event(F x, F y, F dx, F dy) : Event(Point(x, y), Vector(dx, dy)) {}

  /// Evaluates line equation side on given point,

  /// 0 means points lies on the line, -1 left, 1 right
  _ F distance_from(const Point& p) const {
    return p.as_vector().dot(normal) - distance_from_origin;
  }

  /// \brief Return perpendicular event line.
  /// \returns perpendicular event line
  _ Event perpendicular() const {
    return Event(
        this->origin.x, this->origin.y, -this->direction.y, this->direction.x);
  }

  /// Make event translated with given vector.
  _ Event operator+(const Vector& p) const {
    return Event(this->origin + p, this->direction);
  }

  /// Make event translated with given vector.
  _ Event operator-(const Vector& p) const {
    return Event(this->origin - p, this->direction);
  }

  /// Line equation a b coefficients: \f$ a x + b y == c \f$
  const Vector normal;           ///< line equation coefficients \f$ a, b \f$
  const F distance_from_origin;  ///< line equation coefficient \f$ c \f$

  const F b2c;    ///< precalculated \f$ b^2 c \f$
  const F ac;     ///< precalculated \f$ a c \f$
  const F c2;     ///< precalculated \f$ c^2 \f$
  const F inv_b;  ///< precalculated \f$ \frac{1}{b} \f$
};

}  // Barrel
}  // PET2D

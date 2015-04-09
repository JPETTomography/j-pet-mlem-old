#pragma once

#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// Model for 2D barrel PET event

/// Event is described generally by point \f$ (x, y) \f$ and \f$ \phi \f$ angle,
/// however for purpose of various intersection calculations this class holds
/// general line equation \f$ a x + b y + c = 1 \f$ coefficients.
/// This also stores several precalculated variables.
///
/// \note Since \f$ a = sin(\phi), b = -cos(\phi) \f$ then
/// \f$ a^2 + b^2 = 1 \f$.
template <typename FType = double> struct Event : public PET2D::Point<FType> {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Base = Point;

  /// Event requires usage of a concrete constructor.
  Event() = delete;

  /// Make emission event at \f$ (x, y) \f$ point and \f$ \phi \f$ angle.
  _ Event(F x, F y, F phi)
      : Event(x, y, Vector(std::cos(phi), std::sin(phi))) {}

  _ Event(F x, F y, F dx, F dy) : Event(x, y, Vector(dx, dy)) {}

  _ Event(F x, F y, const Vector& direction)
      : Base(x, y),
        direction(direction),
        a(direction.y),
        b(-direction.x),
        // line equation c coefficient: a x + b y == c
        c(a * x + b * y),
        // helper variables
        b2c(b * b * c),
        ac(a * c),
        c2(c * c),
        inv_b(1 / b) {}

 private:
 public:
  /// Make emission event at \f$ p = (x, y) \f$ point and \f$ \phi \f$ angle.
  _ Event(Base p, F phi) : Event(p.x, p.y, phi) {}

  // evaluates line equation side on given point
  // 0 means points lies on the line, -1 left, 1 right
  _ F operator()(const Point& p) const { return a * p.x + b * p.y - c; }

  /// \brief Return perpendicular event line.
  /// \returns perpendicular event line
  _ Event perpendicular() const {
    return Event(this->x, this->y, -direction.y, direction.x);
  }

  /// Make event translated with given vector.
  _ Event operator+(const Vector& p) const {
    return Event(this->x + p.x, this->y + p.y, direction);
  }

  /// Make event translated with given vector.
  _ Event operator-(const Vector& p) const {
    return Event(this->x - p.x, this->y - p.y, direction);
  }

  const Vector direction;

  // const F phi;  ///< \f$ \phi \f$ angle
  /// line equation a b coefficients: a x + b y == c
  const F a;  ///< line equation coefficient \f$ a \f$
  const F b;  ///< line equation coefficient \f$ b \f$
  const F c;  ///< line equation coefficient \f$ c \f$

  const F b2c;    ///< precalculated \f$ b^2 c \f$
  const F ac;     ///< precalculated \f$ a c \f$
  const F c2;     ///< precalculated \f$ c^2 \f$
  const F inv_b;  ///< precalculated \f$ \frac{1}{b} \f$
};
}  // Barrel
}  // PET2D

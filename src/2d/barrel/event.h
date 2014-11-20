#pragma once

#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// Model for 2D barrel PET event
template <typename FType = double> struct Event : public PET2D::Point<FType> {
  using F = FType;
  using Point = PET2D::Point<F>;
  using Base = Point;

  Event() = delete;

  _ Event(F x, F y, F phi) : Base(x, y), phi(phi) {
    // get line equation coefficients
    // a x + b y == c
    a = std::sin(phi);
    b = -std::cos(phi);
    precalculate();
  }

 private:
  _ Event(F x, F y, F phi, F a, F b) : Base(x, y), phi(phi), a(a), b(b) {
    precalculate();
  }

  _ void precalculate() {
    // get line equation coefficients (cont.)
    // a x + b y == c
    c = a * this->x + b * this->y;

    // helper variables
    b2 = b * b;
    b2c = b2 * c;
    ac = a * c;
    a2_b2 = a * a + b2;
    b_a2_b2 = b * a2_b2;
    c2 = c * c;
  }

 public:
  _ Event(Base p, F phi) : Point(p.x, p.y), phi(phi) {}

  // evaluates line equation side on given point
  // 0 means points lies on the line, -1 left, 1 right
  _ F operator()(const Point& p) { return a * p.x + b * p.y - c; }

  _ Event operator+(const Point& p) const {
    return Event(this->x + p.x, this->y + p.y, phi, a, b);
  }

  _ Event operator-(const Point& p) const {
    return Event(this->x - p.x, this->y - p.y, phi, a, b);
  }

  F phi;

  // line equation coefficients
  F a, b, c;
  // precalculated variables
  F b2, b2c, ac, a2_b2, b_a2_b2, c2;
};
}  // Barrel
}  // PET2D

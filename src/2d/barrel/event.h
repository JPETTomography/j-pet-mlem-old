#pragma once

#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// Model for 2D barrel PET event
template <typename FType = double> struct Event {
  using F = FType;
  using Point = PET2D::Point<F>;

  Event(F x, F y, F phi) : x(x), y(y), phi(phi) {
    // get line equation coefficients
    // a x + b y == c
    a = std::sin(phi);
    b = -std::cos(phi);
    precalculate();
  }

 private:
  Event(F x, F y, F phi, F a, F b) : x(x), y(y), phi(phi), a(a), b(b) {
    precalculate();
  }

  void precalculate() {
    // get line equation coefficients (cont.)
    // a x + b y == c
    c = a * x + b * y;

    // helper variables
    b2 = b * b;
    b2c = b2 * c;
    ac = a * c;
    a2_b2 = a * a + b2;
    b_a2_b2 = b * a2_b2;
  }

 public:
  Event(Point p, F phi) : Event(p.x, p.y, phi) {}

  // evaluates line equation side on given point
  // 0 means points lies on the line, -1 left, 1 right
  F operator()(const Point& p) { return a * p.x + b * p.y - c; }

  Event operator+(const Point& p) const {
    return Event(x + p.x, y + p.y, phi, a, b);
  }

  Event operator-(const Point& p) const {
    return Event(x - p.x, y - p.y, phi, a, b);
  }

  F x, y;
  F phi;

  // line equation coefficients
  F a, b, c;
  // precalculated variables
  F b2, b2c, ac, a2_b2, b_a2_b2;
};
}  // Barrel
}  // PET2D

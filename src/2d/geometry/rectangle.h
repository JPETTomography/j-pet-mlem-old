#ifndef RECTANGLE
#define RECTANGLE

#include <iostream>

#include "point.h"
#include "util/random.h"
#include "util/cuda/compat.h"

namespace PET2D {

/** Axis alligned Rectangle centered on center with width = 2*a and height = 2*b
 *
 */
template <typename FType> class Rectangle {
 public:
  using F = FType;
  using Point = PET2D::Point<F>;

  Rectangle(F x, F y, F a, F b) : center(x, y), a(a), b(b), area(4 * a * b) {}

  bool contains(Point p) {
    auto r = p - center;
    return std::abs(r.x) <= a && std::abs(r.y) <= b;
  }

  const Point center;
  const F a;
  const F b;
  const F area;
};
}

#endif  // RECTANGLE

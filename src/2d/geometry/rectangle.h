#pragma once

#include <random>

#include "point.h"

#include "util/cuda/compat.h"

namespace PET2D {

/// Axis alligned & centered rectangle

/// Rectangle that is axis aligned and centered on \c center with $width = 2*a$
/// and $height = 2*b$.
template <typename FType> class Rectangle {
 public:
  using F = FType;
  using Point = PET2D::Point<F>;

  Rectangle(F x, F y, F a, F b) : center(x, y), a(a), b(b), area(4 * a * b) {}

  bool contains(Point p) const {
    auto r = p - center;
    return std::abs(r.x) <= a && std::abs(r.y) <= b;
  }

  const Point center;
  const F a;
  const F b;
  const F area;
};

/// Uniform rectangle point generator
template <typename FType> class RectanglePointGenerator {
 public:
  using F = FType;
  using Rectangle = PET2D::Rectangle<F>;
  using Point = PET2D::Point<F>;

  RectanglePointGenerator(const Rectangle& rec) : rectangle_(rec), uni(-1, 1) {}

  template <class RNG> Point operator()(RNG& rng) {
    F x = uni(rng) * rectangle_.a;
    F y = uni(rng) * rectangle_.b;

    return Point(rectangle_.center.x + x, rectangle_.center.y + y);
  }

 private:
  Rectangle rectangle_;
  std::uniform_real_distribution<F> uni;
};

}  // PET2D

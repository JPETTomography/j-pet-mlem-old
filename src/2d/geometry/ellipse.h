#pragma once

#if !__CUDACC__
#include <iostream>
#endif

#include "point.h"
#include "util/random.h"
#include "util/cuda/compat.h"

namespace PET2D {

/// Unanchored ellipse with given point and axes
template <typename FType> struct Ellipse {
  using F = FType;

  using Point = PET2D::Point<F>;

  Ellipse(F x, F y, F a, F b, F angle)
      : Ellipse(x, y, a, b, angle, compat::sin(angle), compat::cos(angle)) {}

  Ellipse(const Point center, F a, F b, F angle)
      : Ellipse(center.x, center.y, a, b, angle) {}

#if !__CUDACC__
  /// constructs Ellipse from stream
  Ellipse(std::istream& in) : Ellipse(util::read<F>(in), in) {}
#endif

  /// checks if ellipse contains given point
  bool contains(Point p) const {
    auto r = p - center;
    return A * r.x * r.x + 2 * C * r.x * r.y + B * r.y * r.y <= 1;
  }

  const Point center;  ///< ellipse center
  const F a, b;        ///< ellipse axis
  const F angle;       ///< ellipse angle (rotation)
  const F A, B, C;     ///< ellipse equation components
  const F area;        ///< ellipse area

 private:
  Ellipse(F x, F y, F a, F b, F angle, F s, F c)
      : center(x, y),
        a(a),
        b(b),
        angle(angle),
        A(c * c / (a * a) + s * s / (b * b)),
        B(s * s / (a * a) + c * c / (b * b)),
        C(s * c * (1 / (a * a) - 1 / (b * b))),
        area(M_PI * a * b) {}
#if !__CUDACC__
  // Ensure that ellipse params are read in correct order
  Ellipse(F x, std::istream& in) : Ellipse(x, util::read<F>(in), in) {}
  Ellipse(F x, F y, std::istream& in) : Ellipse(x, y, util::read<F>(in), in) {}
  Ellipse(F x, F y, F a, std::istream& in)
      : Ellipse(x, y, a, util::read<F>(in), in) {}
  Ellipse(F x, F y, F a, F b, std::istream& in)
      : Ellipse(x, y, a, b, util::read<F>(in)) {}

 public:
  friend std::ostream& operator<<(std::ostream& out, const Ellipse& ellipse) {
    out << ellipse.center << ' ' << ellipse.a << ' ' << ellipse.b << ' '
        << ellipse.angle;
    return out;
  }
#endif
};

/// Generates random points from given ellipse
template <typename FType> class EllipsePointGenerator {
 public:
  using F = FType;
  using Ellipse = PET2D::Ellipse<F>;
  using Point = PET2D::Point<F>;
  using Distribution = util::random::uniform_real_distribution<F>;

  EllipsePointGenerator(const Ellipse& ellipse)
      : ellipse(ellipse),
        sin(compat::sin(ellipse.angle)),
        cos(compat::cos(ellipse.angle)) {}

  template <class RNG> Point operator()(RNG& rng) {
    F angle = 2 * M_PI * distribution(rng);
    F r = compat::sqrt(distribution(rng));
    F x = ellipse.a * r * std::cos(angle);
    F y = ellipse.b * r * std::sin(angle);

    return Point(cos * x - sin * y + ellipse.center.x,
                 sin * x + cos * y + ellipse.center.y);
  }

 private:
  Ellipse ellipse;
  F sin;
  F cos;
  Distribution distribution;
};

}  // PET2D

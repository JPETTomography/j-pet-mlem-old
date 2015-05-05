#pragma once

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
  Ellipse(std::istream& in)
      : Ellipse(util::read<F>(in),
                util::read<F>(in),
                util::read<F>(in),
                util::read<F>(in),
                util::read<F>(in)) {}
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
};

/// Elliptical emmission source
template <typename FType> struct EllipticalSource : public Ellipse<FType> {
  using F = FType;

  using Ellipse = PET2D::Ellipse<F>;
  using Point = PET2D::Point<F>;

  const F intensity;

  EllipticalSource(F x, F y, F a, F b, F phi, F intensity)
      : Ellipse::Ellipse(x, y, a, b, phi), intensity(intensity) {}

  EllipticalSource(Point center, F a, F b, F phi, F intensity)
      : Ellipse::Ellipse(center, a, b, phi), intensity(intensity) {}

#if !__CUDACC__
  EllipticalSource(std::istream& in)
      : Ellipse::Ellipse(in), intensity(util::read<F>(in)) {}
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

  template <typename Generator> Point operator()(Generator& generator) {
    F angle = 2 * M_PI * distribution(generator);
    F r = compat::sqrt(distribution(generator));
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

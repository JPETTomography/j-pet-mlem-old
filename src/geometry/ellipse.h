#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <random>

#include "point.h"
#include "util/cuda/compat.h"

template <typename FType = double> struct Ellipse {
  typedef FType F;

  Ellipse(F x, F y, F a, F b, F angle)
      : Ellipse(x, y, a, b, angle, compat::sin(angle), compat::cos(angle)) {}

  bool contains(const Point<F> p) const {
    F x = p.x - this->x;
    F y = p.y - this->y;
    return (A * x * x + 2 * C * x * y + B * y * y) <= 1;
  }

  const F x, y;  // center
  const F a, b;
  const F angle;
  const F A, B, C;
  const F measure;

 private:
  Ellipse(F x, F y, F a, F b, F angle, F s, F c)
      : x(x),
        y(y),
        a(a),
        b(b),
        angle(angle),
        A(c * c / (a * a) + s * s / (b * b)),
        B(s * s / (a * a) + c * c / (b * b)),
        C(s * c * (1 / (a * a) - 1 / (b * b))),
        measure(M_PI * a * b) {}
};

template <typename F> class EllipsePointsGenerator {
 public:
  EllipsePointsGenerator(const Ellipse<F>& ellipse)
      : ellipse(ellipse),
        s(compat::sin(ellipse.angle)),
        c(compat::cos(ellipse.angle)) {}

  template <typename G> Point<F> point(G& gen) {
    F angle = 2 * M_PI * uni(gen);
    F r = std::sqrt(uni(gen));
    F x = ellipse.a * r * std::cos(angle);
    F y = ellipse.b * r * std::sin(angle);

    return Point<F>(c * x - s * y + ellipse.x, s * x + c * y + ellipse.y);
  }

 private:
  Ellipse<F> ellipse;
  F s;
  F c;
  std::uniform_real_distribution<F> uni;
};

#endif  // ELLIPSE_H

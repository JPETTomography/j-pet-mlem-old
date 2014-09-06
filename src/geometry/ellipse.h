#ifndef ELLIPSE_H
#define ELLIPSE_H

#include <random>

#include "geometry/point.h"

template <typename F = double> class Ellipse {

 public:
  Ellipse(F x, F y, F a, F b, F angle_rad)
      : x_(x), y_(y), a_(a), b_(b), angle_(angle_rad) {
    F c = std::cos(angle_);
    F s = std::sin(angle_);
    A_ = c * c / (a_ * a_) + s * s / (b_ * b_);
    B_ = s * s / (a_ * a_) + c * c / (b_ * b_);
    C_ = s * c * (F(1.0) / (a_ * a_) - F(1.0) / (b_ * b_));
    measure_ = M_PI * a_ * b_;
  }

  bool in(F x_a, F y_a) {
    F x = x_a - x_;
    F y = y_a - y_;
    return (A_ * x * x + 2 * C_ * x * y + B_ * y * y) <= F(1.0);
  };

  F x() const { return x_; }
  F y() const { return y_; }
  F a() const { return a_; }
  F b() const { return b_; }
  F angle() const { return angle_; }
  F A() const { return A_; }
  F B() const { return B_; }
  F C() const { return C_; }
  F measure() { return measure_; }

 private:
  F x_, y_;  // center
  F a_, b_;
  F angle_;
  F A_, B_, C_;
  F measure_;
};

template <typename F> class EllipsePointsGenerator {
 public:
  EllipsePointsGenerator(const Ellipse<F>& ellipse) : ellipse_(ellipse) {
    c_ = std::cos(ellipse_.angle());
    s_ = std::sin(ellipse_.angle());
  }

  template <typename G> Point<F> point(G& gen) {
    F angle = 2 * M_PI * uni_(gen);
    F r = std::sqrt(uni_(gen));
    F x = ellipse_.a() * r * std::cos(angle);
    F y = ellipse_.b() * r * std::sin(angle);

    return Point<F>(c_ * x - s_ * y + ellipse_.x(),
                    s_ * x + c_ * y + ellipse_.y());
  }

 private:
  Ellipse<F> ellipse_;
  F s_;
  F c_;
  std::uniform_real_distribution<F> uni_;
};

#endif  // ELLIPSE_H

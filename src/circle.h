#pragma once

#include "point.h"
#include "event.h"
#include "svg_ostream.h"

// produces secant angles circle/line intersection as a equation system solution
// see /math/secant.nb
template <typename F = double> class Circle {
 public:
  Circle(F radius)
      : radius_(radius),          // store radius
        radius2_(radius* radius)  // store precomputed square
        {
  }

  typedef F Angle;
  typedef Point<F> Point;
  typedef Event<F> Event;
  typedef std::pair<Point, Point> Secant;

  Secant secant(Event& e) {
    auto sq = sqrt(e.b2 * (-(e.c * e.c) + e.a2_b2 * radius2_));
    auto asq = e.a * sq;

    return std::make_pair(
        Point((e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2),
        Point((e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2));
  }

  F angle(Point p) { return std::atan2(p.y, p.x); }

  std::pair<Angle, Angle> secant_angles(Event& e) {
    auto s = secant(e);
    return std::make_pair(angle(s.first), angle(s.second));
  }

  int section(F angle, int n_detectors) {
    // converting angles to [0,2 Pi) interval
    F normalised_angle = angle > 0 ? angle : (F) 2.0 * M_PI + angle;
    return static_cast<int>(
        round(normalised_angle * n_detectors * INV_TWO_PI)) % n_detectors;
  }

  std::pair<int, int> secant_sections(Event& e, size_t n_detectors) {
    auto sa = secant_angles(e);

    return std::make_pair(section(sa.first, n_detectors),
                          section(sa.second, n_detectors));
  }

  F radius() const { return radius_; }
  F radius2() const { return radius2_; }

  friend svg_ostream<F>& operator<<(svg_ostream<F>& svg, Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius_ << "\"/>" << std::endl;
    return svg;
  }

 private:
  const F radius_;
  const F radius2_;
  static const F TWO_PI;
  static const F INV_TWO_PI;
};

template <typename F> const F Circle<F>::TWO_PI = (F) 2.0 * M_PI;
template <typename F> const F Circle<F>::INV_TWO_PI = (F) 1.0 / TWO_PI;

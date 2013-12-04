#pragma once

#include <initializer_list>

#include "point.h"
#include "2d_xy/event.h"
#include "util/svg_ostream.h"

/// Produces secant angles circle/line intersection as a equation system
/// solution.
/// @see /math/secant.nb
/// @note This circle has only radius specified and center point lies in (0, 0).
template <typename F = double, typename S = int> class Circle {
 public:
  Circle(F radius)
      : radius_(radius),           // store radius
        radius2_(radius * radius)  // store precomputed square
  {}

  typedef F Angle;
  typedef ::Point<F> Point;
  typedef ::Event<F> Event;
  typedef std::initializer_list<Point> Secant;
  typedef std::initializer_list<Angle> SecantAngle;
  typedef std::initializer_list<S> SecantSections;

  Secant secant(const Event& e) {
    auto cabr2 = (-(e.c * e.c) + e.a2_b2 * radius2_);
    auto sq = sqrt(e.b2 * cabr2);
    auto asq = e.a * sq;

    return Secant({ Point((e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2),
                    Point((e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2) });
  }

  F angle(Point p) { return std::atan2(p.y, p.x); }

  SecantAngle secant_angles(Event& e) {
    auto s = secant(e);
    return SecantAngle(angle(s.first), angle(s.second));
  }

  S section(F angle, S n_detectors) {
    // converting angles to [0,2 Pi) interval
    F normalised_angle = angle > 0 ? angle : (F)2.0 * M_PI + angle;
    return static_cast<S>(round(normalised_angle * n_detectors * INV_TWO_PI)) %
           n_detectors;
  }

  SecantSections secant_sections(Event& e, S n_detectors) {
    auto sa = secant_angles(e);

    return SecantSections(section(sa.first, n_detectors),
                          section(sa.second, n_detectors));
  }

  F radius() const { return radius_; }
  F radius2() const { return radius2_; }

  friend svg_ostream<F>& operator<<(svg_ostream<F>& svg, Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius_ << "\"/>" << std::endl;
    return svg;
  }

 private:
  F radius_;
  F radius2_;
  static const F TWO_PI;
  static const F INV_TWO_PI;
};

template <typename F, typename S> const F Circle<F, S>::TWO_PI = (F)2.0 * M_PI;
template <typename F, typename S>
const F Circle<F, S>::INV_TWO_PI = (F)1.0 / TWO_PI;

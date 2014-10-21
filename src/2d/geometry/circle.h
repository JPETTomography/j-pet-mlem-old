#pragma once

#include <initializer_list>

#include "point.h"
#include "2d/barrel/event.h"
#include "util/svg_ostream.h"
#include "util/array.h"

namespace PET2D {

/// Circle with given radius lying in (0, 0)

/// Produces secant angles circle/line intersection as a equation system
/// solution.
///
/// \see
///   /math/secant.nb
/// \note
///   This circle has only radius specified and center point lies in (0, 0).
template <typename FType = double, typename SType = int> class Circle {
 public:
  typedef FType F;
  typedef SType S;

  Circle(F radius)
      : radius_(radius),           // store radius
        radius2_(radius * radius)  // store precomputed square
  {}

  typedef F Angle;
  typedef PET2D::Point<F> Point;
  typedef Barrel::Event<F> Event;
  typedef util::array<2, Point> Secant;
  typedef util::array<2, Angle> SecantAngles;
  typedef util::array<2, S> SecantSections;

  Secant secant(const Event& e) {
    auto cabr2 = (-(e.c * e.c) + e.a2_b2 * radius2_);
    auto sq = sqrt(e.b2 * cabr2);
    auto asq = e.a * sq;

    return Secant({ Point((e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2),
                    Point((e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2) });
  }

  F angle(Point p) { return std::atan2(p.y, p.x); }

  SecantAngles secant_angles(Event& e) {
    SecantAngles sa;
    for (auto& p : secant(e)) {
      sa.push_back(angle(p));
    }
    return sa;
  }

  S section(F angle, S n_detectors) {
    const F TWO_PI = F(2 * M_PI);
    const F INV_TWO_PI = 1 / TWO_PI;

    // converting angles to [0,2 Pi) interval
    F normalised_angle = angle > 0 ? angle : TWO_PI + angle;
    return static_cast<S>(
               std::round(normalised_angle * n_detectors * INV_TWO_PI)) %
           n_detectors;
  }

  SecantSections secant_sections(Event& e, S n_detectors) {
    SecantSections ss;
    for (auto& sa : secant_angles(e)) {
      ss.push_back(section(sa, n_detectors));
    }
    return ss;
  }

  F radius() const { return radius_; }
  F radius2() const { return radius2_; }

  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius_ << "\"/>" << std::endl;
    return svg;
  }

 private:
  F radius_;
  F radius2_;
};
}  // PET2D

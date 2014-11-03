#pragma once

#include <initializer_list>

#include "point.h"
#include "2d/barrel/event.h"
#include "util/svg_ostream.h"
#include "util/array.h"

namespace PET2D {

/// Zero (0, 0) anchored circle with given radius

/// Produces secant angles circle/line intersection as a equation system
/// solution.
///
/// \see
///   /math/secant.nb
/// \note
///   This circle has only radius specified and center point lies in (0, 0).
template <typename FType = double, typename SType = int> class Circle {
 public:
  using F = FType;
  using S = SType;

  Circle(F radius)
      : radius(radius),           // store radius
        radius2(radius * radius)  // store precomputed square
  {}

  Circle& operator=(const Circle& other) { return *(new (this) Circle(other)); }

  using Angle = F;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;
  using Secant = util::array<2, Point>;
  using SecantAngles = util::array<2, Angle>;
  using SecantSections = util::array<2, S>;

  Secant secant(const Event& e) {
    auto cabr2 = (-(e.c * e.c) + e.a2_b2 * radius2);
    auto sq2 = e.b2 * cabr2;
    if (sq2 > 0) {
      auto sq = sqrt(sq2);
      auto asq = e.a * sq;
      return Secant(
          { Point((e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2),
            Point((e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2) });
    } else if (sq2 == 0) {
      return Secant({ Point(e.ac / e.a2_b2, e.b2c / e.b_a2_b2) });
    } else {
      return Secant();
    }
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

  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius << "\"/>" << std::endl;
    return svg;
  }

  const F radius;
  const F radius2;
};
}  // PET2D

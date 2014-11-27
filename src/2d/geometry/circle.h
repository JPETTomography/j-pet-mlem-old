#pragma once

#include "point.h"
#include "2d/barrel/event.h"
#include "util/array.h"
#include "util/cuda/compat.h"
#if !__CUDACC__
#include "util/svg_ostream.h"
#endif

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

  _ Circle(F radius)
      : radius(radius),           // store radius
        radius2(radius * radius)  // store precomputed square
  {}

  // allows copying whole object
  _ Circle& operator=(const Circle& other) { return *new (this) Circle(other); }

  const F radius;
  const F radius2;

  using Angle = F;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;
  using Secant = util::array<2, Point>;
  using SecantAngles = util::array<2, Angle>;
  using SecantSections = util::array<2, S>;

  _ bool intersects(const Event& e) {
    return e.a2_b2 * radius2 > e.c2;
  }

  _ Secant secant(const Event& e) {
    auto cabr2 = e.a2_b2 * radius2 - e.c2;
    if (cabr2 > 0) {
      auto sq = compat::sqrt(e.b2 * cabr2);
      auto asq = e.a * sq;
      return Secant{ Point((e.ac - sq) / e.a2_b2, (e.b2c + asq) / e.b_a2_b2),
                     Point((e.ac + sq) / e.a2_b2, (e.b2c - asq) / e.b_a2_b2) };
    } else if (cabr2 == 0) {
      return Secant{ Point(e.ac / e.a2_b2, e.b2c / e.b_a2_b2) };
    } else {
      return Secant();
    }
  }

  _ F angle(Point p) { return compat::atan2(p.y, p.x); }

  SecantAngles secant_angles(Event& e) {
    SecantAngles sa;
    for (auto& p : secant(e)) {
      sa.push_back(angle(p));
    }
    return sa;
  }

  _ S section(F angle, S n_detectors) {
    const F TWO_PI = F(2 * M_PI);
    const F INV_TWO_PI = 1 / TWO_PI;

    // converting angles to [0,2 Pi) interval
    F normalised_angle = angle > 0 ? angle : TWO_PI + angle;
    return static_cast<S>(
               compat::round(normalised_angle * n_detectors * INV_TWO_PI)) %
           n_detectors;
  }

  SecantSections secant_sections(Event& e, S n_detectors) {
    SecantSections ss;
    for (auto& sa : secant_angles(e)) {
      ss.push_back(section(sa, n_detectors));
    }
    return ss;
  }

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius << "\"/>" << std::endl;
    return svg;
  }
#endif
};
}  // PET2D

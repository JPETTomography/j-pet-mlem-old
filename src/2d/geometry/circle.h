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
template <typename FType> class Circle {
 public:
  using F = FType;

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

  _ bool intersects(const Event& e) const { return radius2 > e.c2; }

  /// \brief Returns secant for given event line
  ///
  /// Event line is described with line equation:
  /// \f$ ax + by + c = 0 \f$ with assumption \f$ a^2 + b^2 = 1 \f$ which is
  /// satisfied for Event since \f$ a = sin(\phi), b = -cos(\phi) \f$.
  _ Secant secant(const Event& e) const {
    auto diff = radius2 - e.c2;
    if (diff > 0) {
      auto sq = e.b * compat::sqrt(diff);
      auto asq = e.a * sq;
      return Secant{ Point(e.ac - sq, (e.b2c + asq) * e.inv_b),
                     Point(e.ac + sq, (e.b2c - asq) * e.inv_b) };
    } else if (diff == 0) {
      return Secant{ Point(e.ac, e.b2c * e.inv_b) };
    } else {
      return Secant();
    }
  }

  _ F angle(Point p) const { return compat::atan2(p.y, p.x); }

  SecantAngles secant_angles(Event& e) const {
    SecantAngles sa;
    for (auto& p : secant(e)) {
      sa.push_back(angle(p));
    }
    return sa;
  }

  template <typename S> _ S section(F angle, S n_detectors) const {
    const F TWO_PI = F(2 * M_PI);
    const F INV_TWO_PI = 1 / TWO_PI;

    // converting angles to [0,2 Pi) interval
    F normalised_angle = angle > 0 ? angle : TWO_PI + angle;
    return static_cast<S>(
               compat::round(normalised_angle * n_detectors * INV_TWO_PI)) %
           n_detectors;
  }

//  SecantSections secant_sections(Event& e, S n_detectors) const {
//    SecantSections ss;
//    for (auto& sa : secant_angles(e)) {
//      ss.push_back(section(sa, n_detectors));
//    }
//    return ss;
//  }

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          Circle& c) {
    svg << "<circle cx=\"0\" cy=\"0\" r=\"" << c.radius << "\"/>" << std::endl;
    return svg;
  }
#endif
};

}  // PET2D

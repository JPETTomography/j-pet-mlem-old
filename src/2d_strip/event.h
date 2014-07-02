#pragma once

#include <cmath>

#include "util/compat.h"

template <typename F> struct Event {
  F z_u;
  F z_d;
  F dl;

  Event(F z_u, F z_d, F dl) $ : z_u(z_u), z_d(z_d), dl(dl) {}

  F tan(const F R) const $ { return (z_u - z_d) / (2 * R); }

  F y(const F tan) const $ {
    return -F(0.5) * (dl / compat::sqrt(1 + tan * tan));
  }

  F z(const F y, const F tan) const $ {
    return F(0.5) * (z_u + z_d + (2 * y * tan));
  }
};

template <typename F> struct ImageSpaceEventTan;

template <typename F> struct ImageSpaceEventAngle {
  const F y;
  const F z;
  const F angle;

  ImageSpaceEventAngle(F y, F z, F angle) $ : y(y), z(z), angle(angle) {}

  ImageSpaceEventTan<F> to_tan() const $ {
    return ImageSpaceEventTan<F>(y, z, compat::tan(angle));
  }
};

template <typename F> struct ImageSpaceEventTan {
  const F y;
  const F z;
  const F tan;

  ImageSpaceEventTan(F y, F z, F tan) : y(y), z(z), tan(tan) {}

  ImageSpaceEventAngle<F> to_angle() const $ {
    return ImageSpaceEventAngle<F>(y, z, compat::atan(tan));
  }
};

template <typename FType> struct EllipseParameters {
  typedef FType F;
  F x, y, a, b;
  F angle;
  F emissions;

  EllipseParameters(F x, F y, F a, F b, F angle, F emissions) $
      : x(x),
        y(y),
        a(a),
        b(b),
        angle(angle),
        emissions(emissions) {}
};

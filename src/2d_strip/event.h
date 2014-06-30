#pragma once

#include <cmath>

template <typename F> struct Event {
  F z_u;
  F z_d;
  F dl;

  Event(F z_u, F z_d, F dl) : z_u(z_u), z_d(z_d), dl(dl) {}

  F tan(const F R) const { return (z_u - z_d) / (2 * R); }

  F y(const F tan) const { return -F(0.5) * (dl / std::sqrt(1 + tan * tan)); }

  F z(const F y, const F tan) const {
    return F(0.5) * (z_u + z_d + (2 * y * tan));
  }
};

template <typename F> struct ImageSpaceEventTan;

template <typename F> struct ImageSpaceEventAngle {
  const F y;
  const F z;
  const F angle;

  ImageSpaceEventAngle(F y, F z, F angle) : y(y), z(z), angle(angle) {}

  ImageSpaceEventTan<F> to_tan() const {
    return ImageSpaceEventTan<F>(y, z, tan(angle));
  }
};

template <typename F> struct ImageSpaceEventTan {
  const F y;
  const F z;
  const F tan;

  ImageSpaceEventTan(F y, F z, F tan) : y(y), z(z), tan(tan) {}

  ImageSpaceEventAngle<F> to_angle() const {
    return ImageSpaceEventAngle<F>(y, z, atan(tan));
  }
};

template <typename T> struct EllipseParameters {
  T x, y, a, b;
  T angle;
  T iter;
};

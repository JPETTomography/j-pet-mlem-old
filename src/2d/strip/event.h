#pragma once

#include <cmath>

#include "util/cuda/compat.h"

namespace PET2D {
namespace Strip {

/// Detector-space event
template <typename FType = double> struct Event {
  using F = FType;
  F z_u;
  F z_d;
  F dl;

  _ Event(F z_u, F z_d, F dl) : z_u(z_u), z_d(z_d), dl(dl) {}

  _ void transform(F R, F& tan, F& y, F& z) const {
    tan = this->tan(R);
    y = this->y(tan);
    z = this->z(y, tan);
  }

 private:
  _ F tan(const F R) const { return (z_u - z_d) / (2 * R); }

  _ F y(const F tan) const {
    return -F(0.5) * (dl / compat::sqrt(1 + tan * tan));
  }

  _ F z(const F y, const F tan) const {
    return F(0.5) * (z_u + z_d + (2 * y * tan));
  }
};

template <typename FType = double> struct ImageSpaceEventTan;

/// Image-space event using angle in radians
template <typename FType = double> struct ImageSpaceEventAngle {
  using F = FType;
  const F y;
  const F z;
  const F angle;

  _ ImageSpaceEventAngle(F y, F z, F angle) : y(y), z(z), angle(angle) {}

  _ ImageSpaceEventTan<F> to_tan() const {
    return ImageSpaceEventTan<F>(y, z, compat::tan(angle));
  }
};

/// Image-space event using angle tangent
template <typename FType> struct ImageSpaceEventTan {
  using F = FType;
  const F y;
  const F z;
  const F tan;

  ImageSpaceEventTan(F y, F z, F tan) : y(y), z(z), tan(tan) {}

  _ ImageSpaceEventAngle<F> to_angle() const {
    return ImageSpaceEventAngle<F>(y, z, compat::atan(tan));
  }
};

/// Describes emissions from ellipse
template <typename FType = double> struct EllipseParameters {
  using F = FType;
  F x, y, a, b;
  F angle;
  F n_emissions;

  _ EllipseParameters(F x, F y, F a, F b, F angle, F n_emissions)
      : x(x), y(y), a(a), b(b), angle(angle), n_emissions(n_emissions) {}
};
}  // Strip
}  // PET2D

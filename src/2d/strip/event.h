#pragma once

#include "util/cuda/compat.h"
#include "2d/geometry/event.h"

namespace PET2D {
namespace Strip {

/// Scanner-space event
template <typename FType> struct Event {
  using F = FType;
  F z_u;
  F z_d;
  F dl;

  _ Event(F z_u, F z_d, F dl) : z_u(z_u), z_d(z_d), dl(dl) {}
  _ Event() {}

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

template <typename FType> struct ImageSpaceEventTan;

/// Image-space event using angle in radians

template <typename FType>
struct ImageSpaceEventAngle : public PET2D::Event<FType> {
  using F = FType;
  using Base = PET2D::Event<F>;
  using Point = PET2D::Point<F>;

  _ ImageSpaceEventAngle(const Base& event) : Base(event) {}
  _ ImageSpaceEventAngle(F y, F z, F angle) : Base(Point(z, y), angle) {}

  _ ImageSpaceEventTan<F> to_tan() const {
    return ImageSpaceEventTan<F>(
        this->center.y, this->center.x, this->direction.x / this->direction.y);
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
template <typename FType> struct EllipseParameters {
  using F = FType;
  F x, y, a, b;
  F angle;
  F n_emissions;

  _ EllipseParameters(F x, F y, F a, F b, F angle, F n_emissions)
      : x(x), y(y), a(a), b(b), angle(angle), n_emissions(n_emissions) {}
};
}  // Strip
}  // PET2D

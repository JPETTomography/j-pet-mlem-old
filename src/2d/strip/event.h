#pragma once

#if !__CUDACC__
#include <ostream>
#endif

#include "util/cuda/compat.h"
#include "2d/geometry/event.h"

namespace PET2D {
namespace Strip {

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
        this->origin.y, this->origin.x, this->direction.x / this->direction.y);
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

}  // Strip
}  // PET2D

#pragma once

#if !__CUDACC__
#include <ostream>
#endif

#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"
#include "2d/barrel/event.h"

#include "util/cuda/compat.h"

namespace PET3D {

/// Generic 3D emission event
////
/// Emission event consist of origin (emission) point \f$ origin = (x, y, z) \f$
/// and direction (vector) \f$ direction = (dx, dy, dz) \f$.
template <typename FType> class Event {
  using Vector = PET3D::Vector<FType>;
  using Vector2D = PET2D::Vector<FType>;
  using Point = PET3D::Point<FType>;
  using BarrelEvent = PET2D::Barrel::Event<FType>;

 public:
  _ Event(const Point& origin, const Vector& direction)
      : origin(origin), direction(direction) {}

  _ BarrelEvent to_barrel_event() const {
    auto direction_2d = Vector2D(direction.x, direction.y);
    direction_2d.normalize();
    return BarrelEvent(origin.x, origin.y, direction_2d.x, direction_2d.y);
  }

  const Point origin;
  const Vector direction;

#if !__CUDACC__
  friend std::ostream& operator<<(std::ostream& out, const Event& event) {
    out << event.origin.x << ' ' << event.origin.y << ' ' << event.origin.z
        << " ";
    out << event.direction.x << ' ' << event.direction.y << " "
        << event.direction.z;
    return out;
  }
#endif
};

}  // PET3D

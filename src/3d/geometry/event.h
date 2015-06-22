#pragma once

#include "3d/geometry/plane.h"
#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"
#include "2d/barrel/event.h"

namespace PET3D {

template <typename FType> class Event {
  using Vector = PET3D::Vector<FType>;
  using Vector2D = PET2D::Vector<FType>;
  using Point = PET3D::Point<FType>;
  using BarrelEvent = PET2D::Barrel::Event<FType>;
  using Plane = PET3D::Plane<FType>;

 public:
  Event(const Point& origin, const Vector& direction)
      : origin(origin), direction(direction) {}

  BarrelEvent to_barrel_event() const {
    auto dir_2d = Vector2D(direction.x, direction.y);
    dir_2d.normalize();
    return BarrelEvent(origin.x, origin.y, dir_2d.x, dir_2d.y);
  }

  Plane z_plane() const { return Plane(1.0, 0.0, 0.0, 1.0); }

  const Point origin;
  const Vector direction;
};

template <typename FType>
std::ostream& operator<<(std::ostream& out, const Event<FType>& ev) {
  auto orig = ev.origin;
  out << orig.x << " " << orig.y << " " << orig.z << " ";
  auto dir = ev.direction;
  out << dir.x << " " << dir.y << " " << dir.z;
  return out;
}
}

#ifndef EVENT_H
#define EVENT_H

#include "geometry/geometry.h"
#include "3d/geometry/plane.h"
#include "2d/barrel/event.h"

namespace PET3D {

template <typename FType> class Event {
  using Vector = Geometry::Vector<3, FType>;
  using Point = Geometry::Point<3, FType>;
  using Event2D = PET2D::Barrel::Event<FType>;
  using Plane = PET3D::Plane<FType>;

 public:
  Event();
  ~Event();

  Event2D to_barrel_event() const {
    return Event2D(
        origin[0], origin[1], std::atan2(direction[1], direction[0]));
  }

  Plane z_plane() const { return Plane(1.0, 0.0, 0.0, 1.0); }

 private:
  const Point origin;
  const Vector direction;
};
}

#endif  // EVENT_H

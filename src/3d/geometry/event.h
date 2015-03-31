#ifndef EVENT_H
#define EVENT_H

#include "geometry/geometry.h"
#include "2d/barrel/event.h"

namespace PET3D {

template <typename FType> class Event {
  using Vector = Geometry::Vector<3, FType>;
  using Point = Geometry::Point<3, FType>;
  using Event2D = PET2D::Barrel::Event<FType>;

 public:
  Event();
  ~Event();

  Event2D to_barrel_event() const {
    FType x = origin[0];
    FType y = origin[1];
    return Event2D(x, y, std::atan2(1.0,0.0));
  }

 private:
  const Point origin;
  const Vector direction;
};
}

#endif  // EVENT_H

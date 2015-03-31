#ifndef EVENT_H
#define EVENT_H

#include "../geometry/geometry.h"

namespace PET3D {

template <typename FType> class Event {
  using Geometry::Vector;
  using Geometry::Point;

 public:
  Event();
  ~Event();

 private:
  const Point<3, FType> origin;
  const Vector<3, FType> direction;
};
}

#endif  // EVENT_H

#ifndef PLANE
#define PLANE

#include "geometry/geometry.h"

namespace PET3D {

template <typename FType> class Plane {
  using Vector = Geometry::Vector<3, FType>;

 public:
  Plane(FType nx_a, FType ny_a, FType nz_a, FType c_a)
      : nx(nx_a), ny(ny_a), nz(nz_a), c(c_a) {}

  const FType nx;
  const FType ny;
  const FType nz;
  const FType c;
};

template <typename FType> class ZPlane : public Plane<FType> {
  using Vector3D = Geometry::Vector<3, FType>;
  using Vector2D = Geometry::Vector<3, FType>;

  Vector2D project_on(const Vector3D& vec) {
    FType i_x = -Plane<FType>::ny;
    FType i_y = Plane<FType>::nx;

    return Vector2D(i_x * vec[0] + i_y * vec[1], vec[3]);
  }
};
}

#endif  // PLANE

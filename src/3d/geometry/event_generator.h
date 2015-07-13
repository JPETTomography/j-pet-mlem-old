#pragma once

#include "event.h"
#include "vector.h"
#include "point.h"
#include "distribution.h"

namespace PET3D {

/// Generates events from given voxel
template <typename FType> class VoxelEventGenerator {
 public:
  using F = FType;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;

  VoxelEventGenerator(const Point& lover_left_corner, const Vector& size)
      : lover_left_corner(lover_left_corner),
        uni_x(0, size.x),
        uni_y(0, size.y),
        uni_z(0, size.z) {}

  template <class RNG> Event operator()(RNG& rng) {
    F x = lover_left_corner.x + uni_x(rng);
    F y = lover_left_corner.y + uni_y(rng);
    F z = lover_left_corner.z + uni_z(rng);

    return Event(Point(x, y, z), spherical_distribution(rng));
  }

 private:
  const Point lover_left_corner;
  std::uniform_real_distribution<F> uni_x;
  std::uniform_real_distribution<F> uni_y;
  std::uniform_real_distribution<F> uni_z;
  PET3D::Distribution::SphericalDistribution<F> spherical_distribution;
};

}  // PET2D

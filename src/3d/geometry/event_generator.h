#ifndef EVENT_GENERATOR
#define EVENT_GENERATOR

#include <random>

#include "event.h"
#include "vector.h"
#include "point.h"

namespace PET3D {
template <typename F = float> class spherical_distribution {
 public:
  using Vector = PET3D::Vector<F>;

  spherical_distribution(F R = 1)
      : R(R), R2(R * R), phi_dist(-M_PI, M_PI), z_dist(-R, R) {}
  const F R;
  const F R2;

  template <typename RNG> Vector operator()(RNG& rng) {
    F phi = phi_dist(rng);
    F z = z_dist(rng);
    F r = std::sqrt(R2 - z * z);
    F x = r * cos(phi);
    F y = r * sin(phi);
    return Vector(x, y, z);
  }

 private:
  std::uniform_real_distribution<F> phi_dist;
  std::uniform_real_distribution<F> z_dist;
};

template <typename F> class voxel_event_generator {
 public:
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;

  voxel_event_generator(const Point& lover_left_corner, const Vector& size)
      : lover_left_corner(lover_left_corner),
        uni_x(0, size.x),
        uni_y(0, size.y),
        uni_z(0, size.z) {}

  template <typename RNG> Event operator()(RNG& rng) {
    F x = lover_left_corner.x + uni_x(rng);
    F y = lover_left_corner.y + uni_y(rng);
    F z = lover_left_corner.z + uni_z(rng);

    return Event(
                Point(x, y, z),
                spherical_distribution(rng)
                );
  }

 private:
  const Point lover_left_corner;
  std::uniform_real_distribution<F> uni_x;
  std::uniform_real_distribution<F> uni_y;
  std::uniform_real_distribution<F> uni_z;
  PET3D::spherical_distribution<F>  spherical_distribution;
};
}
#endif  // EVENT_GENERATOR

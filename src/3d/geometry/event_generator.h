#ifndef EVENT_GENERATOR
#define EVENT_GENERATOR

#include <random>

#include "event.h"
#include "vector.h"
#include "point.h"

namespace PET3D {
template <typename F> class spherical_distribution {
 public:
  using Vector = PET3D::Vector<F>;

  spherical_distribution(F theta_min = -M_PI / 2, F theta_max = M_PI / 2)
      : phi_dist(-M_PI, M_PI), z_dist(sin(theta_min), sin(theta_max)) {}

  template <typename RNG> Vector operator()(RNG& rng) {

    F z = z_dist(rng);
    F r = std::sqrt(1 - z * z);
    F phi = phi_dist(rng);
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

    return Event(Point(x, y, z), spherical_distribution(rng));
  }

 private:
  const Point lover_left_corner;
  std::uniform_real_distribution<F> uni_x;
  std::uniform_real_distribution<F> uni_y;
  std::uniform_real_distribution<F> uni_z;
  PET3D::spherical_distribution<F> spherical_distribution;
};

template <typename F> class cylinder_point_generator {
 public:
  using Point = PET3D::Point<F>;
  cylinder_point_generator(F radius, F height)
      : radius(radius),
        height(height),
        uni_h(-height / 2, height / 2),
        uni_phi(0, 2 * M_PI),
        uni_r(0, 1) {}

  template <typename RNG> Point operator()(RNG& rng) {
    F phi = uni_phi(rng);
    F r = radius * std::sqrt(uni_r(rng));
    F x = r * cos(phi);
    F y = r * sin(phi);
    F h = uni_h(rng);
    return Point(x, y, h);
  }

 private:
  const F radius;
  const F height;
  std::uniform_real_distribution<F> uni_h;
  std::uniform_real_distribution<F> uni_phi;
  std::uniform_real_distribution<F> uni_r;
};
}
#endif  // EVENT_GENERATOR

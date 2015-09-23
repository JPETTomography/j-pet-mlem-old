#pragma once

#include "event.h"
#include "vector.h"
#include "point.h"

#include "util/cuda/compat.h"
#include "util/random.h"

#if !__CUDACC__
#include "util/json.h"
#endif

namespace PET3D {
/// Various 3D distribution types
namespace Distribution {

template <typename FType> class SphericalDistribution {
 public:
  using F = FType;
  using Vector = PET3D::Vector<F>;

  _ SphericalDistribution(F theta_min = F(-M_PI) / 2, F theta_max = F(M_PI) / 2)
      : theta_min(theta_min),
        theta_max(theta_max),
        phi_dist(F(-M_PI), F(M_PI)),
        z_dist(sin(theta_min), sin(theta_max)) {}

#if !__CUDACC__
  SphericalDistribution(const json& j)
      : SphericalDistribution(
            j.count("theta-min") ? j["theta-min"].get<F>() : F(-M_PI) / 2,
            j.count("theta-max") ? j["theta-max"].get<F>() : F(M_PI) / 2) {
    std::string type = j["type"];
    if (type != "spherical") {
      throw("type mismatch in SphericalDistribution: " + type);
    }
  }
#endif

  template <class RNG> _ Vector operator()(RNG& rng) {
    F z = z_dist(rng);
    F r = compat::sqrt(1 - z * z);
    F phi = phi_dist(rng);
    F x = r * cos(phi);
    F y = r * sin(phi);
    return Vector(x, y, z);
  }

  const F theta_min;
  const F theta_max;

 private:
  util::random::uniform_real_distribution<F> phi_dist;
  util::random::uniform_real_distribution<F> z_dist;
};

template <typename FType> class SingleDirectionDistribution {
 public:
  using F = FType;
  using Vector = PET3D::Vector<F>;
  SingleDirectionDistribution(const Vector& direction)
      : direction(direction.normalized()) {}

#if !__CUDACC__
  SingleDirectionDistribution(const json& j)
      : SingleDirectionDistribution(Vector(j["direction"])) {
    std::string type = j["type"];
    if (type != "single-direction") {
      throw("type mismatch in SingleDirectionDistribution: " + type);
    }
  }
#endif

  template <class RNG> Vector operator()(RNG& rng) {
    (void)rng;  // unused
    return direction;
  }

  const Vector direction;
};

template <typename FType> class CylinderPointDistribution {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  CylinderPointDistribution(F radius, F height)
      : radius(radius),
        height(height),
        uni_h(-height / 2, height / 2),
        uni_phi(0, 2 * F(M_PI)),
        uni_r(0, 1) {}

  template <class RNG> Point operator()(RNG& rng) {
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
  util::random::uniform_real_distribution<F> uni_h;
  util::random::uniform_real_distribution<F> uni_phi;
  util::random::uniform_real_distribution<F> uni_r;
};

template <typename FType> class BallPointDistribution {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  BallPointDistribution(F radius = 1) : radius(radius), uni(-radius, radius) {}

  template <class RNG> Point operator()(RNG& rng) {
    F x, y, z;
    do {
      x = uni(rng);
      y = uni(rng);
      z = uni(rng);
    } while ((x * x + y * y + z * z) > F(1.0));
    return Point(x, y, z);
  }

  const F radius;

 private:
  util::random::uniform_real_distribution<F> uni;
};

template <typename FType> class EllipsoidPointDistribution {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  EllipsoidPointDistribution(F rx, F ry, F rz)
      : rx(rx), ry(ry), rz(rz), ball(1) {}

  const F rx, ry, rz;

  template <class RNG> Point operator()(RNG& rng) {
    Point p = ball(rng);
    p.x *= rx;
    p.y *= ry;
    p.z *= rz;
    return p;
  }

 private:
  BallPointDistribution<F> ball;
};

}  // Distribution
}  // PET2D

#pragma once

#include "3d/geometry/point.h"
#include "3d/geometry/event.h"
#include "3d/geometry/event_generator.h"

template <typename FType, typename RNG, typename AngularDistribution>
class PhantomRegion {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;

  PhantomRegion(F intensity,
                AngularDistribution angular = AngularDistribution())
      : intensity(intensity), angular_distribution(angular){};
  virtual bool in(const Point&) = 0;
  virtual Point random_point(RNG&) = 0;
  Event random_event(RNG& rng) {
    return Event(random_point(rng), angular_distribution(rng));
  }

  Point operator()(RNG& rng) { return random_event(rng); };

  virtual F volume() const = 0;

  const F intensity;

 protected:
  AngularDistribution angular_distribution;
};

template <typename FType,
          typename RNG,
          typename AngularDistribution = PET3D::spherical_distribution<FType>>
class CylinderRegion : public PhantomRegion<FType, RNG, AngularDistribution> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;

  CylinderRegion(F radius, F height, F intensity, AngularDistribution angular)
      : PhantomRegion<FType, RNG, AngularDistribution>(intensity, angular),
        radius(radius),
        height(height),

        dist(radius, height){};

  bool in(const Point& p) {
    return ((p.x * p.x + p.y * p.y < radius * radius) && (p.z <= height / 2) &&
            (p.z >= -height / 2));
  }
  Point random_point(RNG& rng) { return dist(rng); }
  F volume() const { return M_PI * radius * radius * height; }

  const F radius;
  const F height;

 private:
  PET3D::cylinder_point_distribution<F> dist;
};

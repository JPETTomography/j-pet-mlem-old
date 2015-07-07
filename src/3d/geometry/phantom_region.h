#ifndef PHANTOM_REGION
#define PHANTOM_REGION

#include "3d/geometry/point.h"
#include "3d/geometry/event.h"
#include "3d/geometry/event_generator.h"
#include "3d/geometry/matrix.h"

namespace PET3D {

template <typename FType, typename RNG> class PhantomRegion {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;

  PhantomRegion(F intensity) : intensity(intensity) {}
  virtual bool in(const Point&) const = 0;
  virtual Point random_point(RNG&) = 0;
  virtual Vector random_direction(RNG& rng) = 0;

  Event random_event(RNG& rng) {
    return Event(this->random_point(rng), this->random_direction(rng));
  }

  Point operator()(RNG& rng) { return random_event(rng); }

  virtual F volume() const = 0;
  F weight() const { return intensity * volume(); }
  const F intensity;
};

template <typename FType, typename RNG, typename AngularDistribution>
class AbstractPhantomRegion : public PhantomRegion<FType, RNG> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;

  AbstractPhantomRegion(F intensity, AngularDistribution angular)
      : PhantomRegion<F, RNG>(intensity), angular_distribution(angular) {}
  Vector random_direction(RNG& rng) { return angular_distribution(rng); }

 protected:
  AngularDistribution angular_distribution;
};

template <typename FType,
          typename RNG,
          typename AngularDistribution = PET3D::SphericalDistribution<FType>>
class CylinderRegion
    : public AbstractPhantomRegion<FType, RNG, AngularDistribution> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;

  CylinderRegion(
      F radius,
      F height,
      F intensity,
      AngularDistribution angular = PET3D::SphericalDistribution<F>())
      : AbstractPhantomRegion<FType, RNG, AngularDistribution>(intensity,
                                                               angular),
        radius(radius),
        height(height),

        dist(radius, height) {}

  bool in(const Point& p) const {
    return ((p.x * p.x + p.y * p.y < radius * radius) && (p.z <= height / 2) &&
            (p.z >= -height / 2));
  }
  Point random_point(RNG& rng) { return dist(rng); }
  F volume() const { return M_PI * radius * radius * height; }

  const F radius;
  const F height;

 private:
  PET3D::CylinderPointDistribution<F> dist;
};

template <typename FType,
          typename RNG,
          typename AngularDistribution = PET3D::SphericalDistribution<FType>>
class EllipsoidRegion
    : public AbstractPhantomRegion<FType, RNG, AngularDistribution> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;

  EllipsoidRegion(
      F rx,
      F ry,
      F rz,
      F intensity,
      AngularDistribution angular = PET3D::SphericalDistribution<F>())
      : AbstractPhantomRegion<FType, RNG, AngularDistribution>(intensity,
                                                               angular),
        rx(rx),
        ry(ry),
        rz(rz),
        dist(rx, ry, rz) {}

  bool in(const Point& p) const {
    F x = p.x / rx;
    F y = p.y / ry;
    F z = p.z / rz;
    return (x * x + y * y + z * z <= 1.0);
  }

  Point random_point(RNG& rng) { return dist(rng); }
  F volume() const { return 4.0 / 3.0 * M_PI * rx * ry * rz; }

  const F rx, ry, rz;

 private:
  PET3D::EllipsoidPointDistribution<F> dist;
};

/* Rotated -----------------------------------------------------------------
 */

template <typename FType, typename RNG>
class RotatedPhantomRegion : public PhantomRegion<FType, RNG> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;

  RotatedPhantomRegion(PhantomRegion<F, RNG>* region, const Matrix<FType>& R)
      : PhantomRegion<F, RNG>(region->intensity),
        region(region),
        R(R),
        transposed_R(transpose(R)) {}

  F volume() const { return region->volume(); }
  Point random_point(RNG& rng) {
    Point p = region->random_point(rng);
    return Point::from_vector(R * p.as_vector());
  }

  Vector random_direction(RNG& rng) {
    Vector v = region->random_direction(rng);
    return R * v;
  }

  bool in(const Point& p) const {
    return region->in(Point::from_vector(transposed_R * p.as_vector()));
  }

 private:
  PhantomRegion<F, RNG>* region;
  Matrix<F> R;
  Matrix<F> transposed_R;
};

template <typename FType, typename RNG>
class TranslatedPhantomRegion : public PhantomRegion<FType, RNG> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;

  TranslatedPhantomRegion(PhantomRegion<F, RNG>* region,
                          const Vector displacement)
      : PhantomRegion<F, RNG>(region->intensity),
        region(region),
        displacement(displacement) {}

  F volume() const { return region->volume(); }
  Point random_point(RNG& rng) {
    return region->random_point(rng) + displacement;
  }

  Vector random_direction(RNG& rng) { return region->random_direction(rng); }

  bool in(const Point& p) const { return region->in(p - displacement); }

 private:
  PhantomRegion<F, RNG>* region;
  Vector displacement;
};

template <typename FType,
          typename RNG,
          typename AngularDistribution = PET3D::SphericalDistribution<FType>>
class PointRegion
    : public AbstractPhantomRegion<FType, RNG, AngularDistribution> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;
  using Vector = PET3D::Vector<F>;
  PointRegion(F intensity, AngularDistribution angular, const Point& origin)
      : AbstractPhantomRegion<FType, RNG, AngularDistribution>(intensity,
                                                               angular),
        origin(origin) {}

  Point random_point(RNG& rng) {
    (void)rng;  // unused
    return origin;
  }

  bool in(const Point& p) const { return p == origin; }

  F volume() const { return F(1); }

  const Point origin;
};
}

#endif  // PHANTOM_REGION

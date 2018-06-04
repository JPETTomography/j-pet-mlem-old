#pragma once

#include <random>

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

#include "event.h"
#include "event_generator.h"
#include "point.h"
#include "vector.h"
#include "matrix.h"

namespace PET3D {

/// Virtual phantom made of regions
template <class RNGClass, typename FType> class Phantom {
 public:
  using RNG = RNGClass;
  using F = FType;
  using Event = PET3D::Event<F>;
  using Point = PET3D::Point<F>;
  using Vector = PET3D::Vector<F>;

  /// Abstract phantom region (must subclass)
  ////
  /// Must provide at least intensity for the region.
  struct Region {
    Region(F intensity) : intensity(intensity) {}
    virtual bool contains(const Point&) const = 0;
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

  /// Abstract phantom region with angular distribution
  template <typename AngularDistribution>
  class AngularDistributionRegion : public Region {
   public:
    AngularDistributionRegion(F intensity, AngularDistribution angular)
        : Region(intensity), angular_distribution(angular) {}
    Vector random_direction(RNG& rng) { return angular_distribution(rng); }

   protected:
    AngularDistribution angular_distribution;
  };

  /// Cylindrical region
  template <typename AngularDistribution =
                PET3D::Distribution::SphericalDistribution<FType>>
  class CylinderRegion : public AngularDistributionRegion<AngularDistribution> {
   public:
    CylinderRegion(F radius,
                   F height,
                   F intensity,
                   AngularDistribution angular = AngularDistribution())
        : AngularDistributionRegion<AngularDistribution>(intensity, angular),
          radius(radius),
          height(height),
          distribution(radius, height) {}

    bool contains(const Point& p) const {
      return ((p.x * p.x + p.y * p.y < radius * radius) &&
              (p.z <= height / 2) && (p.z >= -height / 2));
    }
    Point random_point(RNG& rng) { return distribution(rng); }
    F volume() const { return M_PI * radius * radius * height; }

    const F radius;
    const F height;

   private:
    PET3D::Distribution::CylinderPointDistribution<F> distribution;
  };

  /// Ellipsoid region
  template <typename AngularDistribution =
                PET3D::Distribution::SphericalDistribution<FType>>
  class EllipsoidRegion
      : public AngularDistributionRegion<AngularDistribution> {
   public:
    EllipsoidRegion(F rx,
                    F ry,
                    F rz,
                    F intensity,
                    AngularDistribution angular = AngularDistribution())
        : AngularDistributionRegion<AngularDistribution>(intensity, angular),
          rx(rx),
          ry(ry),
          rz(rz),
          distribution(rx, ry, rz) {}

    bool contains(const Point& p) const {
      F x = p.x / rx;
      F y = p.y / ry;
      F z = p.z / rz;
      return (x * x + y * y + z * z <= 1.0);
    }

    Point random_point(RNG& rng) { return distribution(rng); }
    F volume() const { return 4.0 / 3.0 * M_PI * rx * ry * rz; }

    const F rx, ry, rz;

   private:
    PET3D::Distribution::EllipsoidPointDistribution<F> distribution;
  };

  template <typename AngularDistribution =
                PET3D::Distribution::SphericalDistribution<FType>>
  class RectangularRegion
      : public AngularDistributionRegion<AngularDistribution> {
   public:
    RectangularRegion(F ax, F ay, F az, F intensity)
        : AngularDistributionRegion<AngularDistribution>(intensity,
                                                         AngularDistribution()),
          ax(ax),
          ay(ay),
          az(az),
          dx(-ax / 2, ax / 2),
          dy(-ay / 2, ay / 2),
          dz(-az / 2, az / 2) {}

    bool contains(const Point& point) const {
      return (point.x <= ax / 2) && (point.x >= -ax / 2) &&
             (point.y <= ay / 2) && point.y >= -ay / 2 && (point.z <= az / 2) &&
             (point.z >= -az / 2);
    }

    Point random_point(RNG& rng) { return Point(dx(rng), dy(rng), dz(rng)); }

    F volume() const { return ax * ay * az; }

    const F ax;
    const F ay;
    const F az;

   private:
    std::uniform_real_distribution<F> dx;
    std::uniform_real_distribution<F> dy;
    std::uniform_real_distribution<F> dz;
  };

  /// Region rotated using given 3x3 transformation matrix
  class RotatedRegion : public Region {
   public:
    using Matrix = PET3D::Matrix<F>;

    RotatedRegion(Region& region, const Matrix& R)
        : Region(region.intensity),
          region(region),
          R(R),
          transposed_R(transpose(R)) {}
    RotatedRegion(Region* region, const Matrix& R)
        : Region(region->intensity),
          region(*region),
          R(R),
          transposed_R(transpose(R)) {}

    F volume() const { return region.volume(); }
    Point random_point(RNG& rng) {
      Point p = region.random_point(rng);
      return Point::from_vector(R * p.as_vector());
    }

    Vector random_direction(RNG& rng) {
      Vector v = region.random_direction(rng);
      return R * v;
    }

    bool contains(const Point& p) const {
      return region.contains(Point::from_vector(transposed_R * p.as_vector()));
    }

   private:
    Region& region;
    const Matrix R;
    const Matrix transposed_R;
  };

  /// Region translated using given vector
  class TranslatedRegion : public Region {
   public:
    TranslatedRegion(Region& region, const Vector displacement)
        : Region(region.intensity),
          region(region),
          displacement(displacement) {}
    TranslatedRegion(Region* region, const Vector displacement)
        : Region(region->intensity),
          region(*region),
          displacement(displacement) {}

    F volume() const { return region.volume(); }
    Point random_point(RNG& rng) {
      return region.random_point(rng) + displacement;
    }

    Vector random_direction(RNG& rng) { return region.random_direction(rng); }

    bool contains(const Point& p) const {
      return region.contains(p - displacement);
    }

   private:
    Region& region;
    const Vector displacement;
  };

  /// Point region
  template <typename AngularDistribution =
                PET3D::Distribution::SphericalDistribution<FType>>
  class PointRegion : public AngularDistributionRegion<AngularDistribution> {
   public:
    using F = FType;
    using Point = PET3D::Point<F>;
    using Event = PET3D::Event<F>;
    using Vector = PET3D::Vector<F>;
    PointRegion(F intensity, AngularDistribution angular, const Point& origin)
        : AngularDistributionRegion<AngularDistribution>(intensity, angular),
          origin(origin) {}

    Point random_point(RNG& rng) {
      (void)rng;  // unused
      return origin;
    }

    bool contains(const Point& p) const { return p == origin; }

    F volume() const { return F(1); }

    const Point origin;
  };

  using RegionPtrList = std::vector<Region*>;

 private:
  RegionPtrList region_list;
  std::vector<F> CDF;

  std::vector<Event> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

  util::random::uniform_real_distribution<F> uniform;
  util::random::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const RegionPtrList& el, bool additive = false)
      : region_list(el),
        CDF(el.size(), 0),
        uniform_angle(-1, 1),
        additive(additive) {
    CDF[0] = region_list[0]->weight();

    for (size_t i = 1; i < el.size(); i++) {
      CDF[i] = region_list[i]->weight() + CDF[i - 1];
    }
    F norm = CDF[el.size() - 1];
    for (size_t i = 0; i < el.size(); i++) {
      CDF[i] /= norm;
#if DEBUG
      std::cerr << "CDF[" << i << "]=" << CDF[i];
#endif
    }
  }

  bool additive;

  size_t n_events() { return events.size(); }

  template <class RNG> size_t choose_region(RNG& rng) {
    F r = uniform(rng);
    size_t i = 0;

    while (r > CDF[i])
      ++i;

    return i;
  }

  Point gen_point(RNG& rng, size_t& region_index) {
  again:
    region_index = choose_region(rng);
    Point p = region_list[region_index]->random_point(rng);
    if (additive)
      return p;
    for (size_t i = 0; i < region_index; i++) {
      if (region_list[i]->contains(p))
        goto again;
    }
    return p;
  }

  Point gen_point(RNG& rng) {
    size_t region_index;
    return gen_point(rng, region_index);
  }

  Event gen_event(RNG& rng) {
    size_t region_index;
    Point p = gen_point(rng, region_index);
    return Event(p, region_list[region_index]->random_direction(rng));
  }
};
}

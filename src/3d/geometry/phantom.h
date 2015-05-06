#pragma once

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

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
          typename AngularDistribution = PET3D::spherical_distribution<FType>>
class CylinderRegion
    : public AbstractPhantomRegion<FType, RNG, AngularDistribution> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;

  CylinderRegion(
      F radius,
      F height,
      F intensity,
      AngularDistribution angular = PET3D::spherical_distribution<F>())
      : AbstractPhantomRegion<FType, RNG, AngularDistribution>(intensity,
                                                               angular),
        radius(radius),
        height(height),

        dist(radius, height){};

  bool in(const Point& p) const {
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

template <typename FType, typename SType, typename RNG> class Phantom {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;

 private:
  std::vector<PhantomRegion<F, RNG>*> region_list;
  std::vector<F> CDF;

  std::vector<Event<F>> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const std::vector<PhantomRegion<F, RNG>*>& el)
      : region_list(el), CDF(el.size(), 0), uniform_angle(-1, 1) {
    CDF[0] = region_list[0]->weight();

    for (size_t i = 1; i < el.size(); i++) {
      CDF[i] = region_list[i]->weight() + CDF[i - 1];
    }
    F norm = CDF[el.size() - 1];
    for (size_t i = 0; i < el.size(); i++) {
      CDF[i] /= norm;
    }
  }

  size_t n_events() { return events.size(); }

  template <typename G> size_t choose_region(G& gen) {
    F r = uniform(gen);
    size_t i = 0;

    while (r > CDF[i])
      ++i;

    return i;
  }

  template <typename Generator> Point<F> gen_point(Generator& generator) {
  again:
    size_t i_region = choose_region(generator);
    Point<F> p = region_list[i_region]->random_point();
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j].shape.contains(p))
        goto again;
    }
    return p;
  }

  template <typename Generator> Event<F> gen_event(Generator& generator) {
  again:
    size_t i_region = choose_region(generator);
    Point<F> p = region_list[i_region]->random_point(generator);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j]->in(p))
        goto again;
    }
    return Event<F>(p, region_list[i_region]->random_direction(generator));
  }
};

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
    return from_vector(R * p.as_vector());
  }

  Vector random_direction(RNG& rng) {
    Vector v = region->random_direction(rng);
    return R * v;
  }

  bool in(const Point& p) const {
    return region->in(from_vector(transposed_R * p.as_vector()));
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
        displacement(displacement),
        region(region) {}

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
}

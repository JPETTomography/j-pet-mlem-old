#pragma once

#include "3d/geometry/point.h"
#include "3d/geometry/event.h"
#include "3d/geometry/event_generator.h"

namespace PET3D {

template <typename FType, typename RNG> class PhantomRegion {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;

  PhantomRegion(F intensity) : intensity(intensity){};
  virtual bool in(const Point&) = 0;
  virtual Point random_point(RNG&) = 0;
  virtual Event random_event(RNG& rng) = 0;

  Point operator()(RNG& rng) { return random_event(rng); };

  virtual F volume() const = 0;

  const F intensity;

 protected:
};

template <typename FType, typename RNG, typename AngularDistribution>
class AbstractPhantomRegion : public PhantomRegion<FType, RNG> {
 public:
  using F = FType;
  using Point = PET3D::Point<F>;
  using Event = PET3D::Event<F>;

  AbstractPhantomRegion(F intensity, AngularDistribution angular)
      : PhantomRegion<F, RNG>(intensity), angular_distribution(angular) {}

  Event random_event(RNG& rng) {
    return Event(this->random_point(rng), angular_distribution(rng));
  }

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

  CylinderRegion(F radius, F height, F intensity, AngularDistribution angular)
      : AbstractPhantomRegion<FType, RNG, AngularDistribution>(intensity,
                                                               angular),
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
#if 0
template <typename FType, typename SType, typename RNG> class Phantom {
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using rng = std::minstd_rand0;

 private:
  std::vector<PhantomRegion<F>> region_list;
  std::vector<F> CDF;
  std::vector<EllipsePointGenerator<F>> point_generators;

  std::vector<Event<F>> events;
  std::vector<std::vector<F>> output;
  std::vector<std::vector<F>> output_without_errors;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const D& detector, const std::vector<PhantomRegion<F>>& el)
      : detector(detector),
        region_list(el),
        CDF(el.size(), 0),
        uniform_angle(-1, 1) {
    CDF[0] = region_list[0].weight;

    for (size_t i = 1; i < el.size(); i++) {
      CDF[i] = region_list[i].weight + CDF[i - 1];
    }
    F norm = CDF[el.size() - 1];
    for (size_t i = 0; i < el.size(); i++) {
      CDF[i] /= norm;
    }

    for (size_t i = 0; i < el.size(); ++i)
      point_generators.emplace_back(el[i].shape);

    output.assign(detector.n_y_pixels, std::vector<F>(detector.n_z_pixels, 0));
    output_without_errors.assign(detector.n_y_pixels,
                                 std::vector<F>(detector.n_z_pixels, 0));
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
    Point<F> p = point_generators[i_region](generator);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j].shape.contains(p))
        goto again;
    }
    return p;
  }

  template <typename Generator>
  ImageSpaceEventAngle<F> gen_event(Generator& generator) {
    Point<F> p = gen_point(generator);
    F rangle = F(M_PI_4) * uniform_angle(generator);
    return ImageSpaceEventAngle<F>(p.y, p.x, rangle);
  }

  void operator()(int n_emissions) {

    std::vector<std::vector<Event<F>>> event_list_per_thread(
        omp_get_max_threads());

    rng rd;
    std::vector<rng> rng_list;

    for (int i = 0; i < omp_get_max_threads(); ++i) {

      rng_list.push_back(rd);
      rng_list[i].seed(42 + (3453 * i));
      // OR
      // Turn on leapfrogging with an offset that depends on the task id
    }

    for (int i = 0; i < omp_get_max_threads(); ++i)
      event_list_per_thread[i].clear();

#if _OPENMP
#pragma omp for schedule(static)
#endif
    for (int emission = 0; emission < n_emissions; ++emission) {

      auto event = gen_event(rng_list[omp_get_thread_num()]);

      auto res = detector.detect_event(event, rng_list[omp_get_thread_num()]);
      if (res.second) {

        ImageSpaceEventTan<F> revent =
            detector.from_projection_space_tan(res.first);

        if (std::abs(revent.y) >= detector.radius)
          continue;

        Pixel pp = detector.pixel_at(Point<F>(event.z, event.y));

        if (detector.contains_pixel(pp)) {

          Pixel p = detector.pixel_at(Point<F>(revent.z, revent.y));

          output[p.y][p.x]++;
          output_without_errors[pp.y][pp.x]++;

          event_list_per_thread[omp_get_thread_num()].push_back(res.first);
        }
      }
    }

    for (int i = 0; i < omp_get_max_threads(); ++i) {
      events.insert(events.end(),
                    event_list_per_thread[i].begin(),
                    event_list_per_thread[i].end());
    }
  }
};
#endif
}

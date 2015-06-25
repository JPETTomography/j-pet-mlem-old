#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>

#include "2d/geometry/event.h"
#include "2d/geometry/ellipse.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

namespace PET2D {
namespace Strip {

// template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

/// Virtual phantom region made of ellipse and intensity
template <typename FType> struct PhantomRegion {
  using F = FType;

  PhantomRegion(const Ellipse<F>& ellipse, F intensity)
      : shape(ellipse), intensity(intensity), weight(intensity * shape.area) {}

  bool contains(Point<F> p) const { return shape.constains(p); }

  const Ellipse<F> shape;
  const F intensity;
  const F weight;
};

/// Virtual phantom made of elliptical regions
template <typename FType, typename SType> class Phantom {
 public:
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using RNG = std::minstd_rand0;
  using Event = PET2D::Event<F>;

 private:
  int n_events_;

  std::vector<PhantomRegion<F>> region_list;
  std::vector<F> CDF;
  std::vector<EllipsePointGenerator<F>> point_generators;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom(const std::vector<PhantomRegion<F>>& el)
      : region_list(el), CDF(el.size(), 0), uniform_angle(-1, 1) {

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
  }

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
  PET2D::Event<F> gen_event(Generator& generator) {
    Point<F> p = gen_point(generator);
    F rangle = F(M_PI_2) * uniform_angle(generator);
    return PET2D::Event<F>(p, rangle);
  }
};
}  // Strip
}  // PET2D

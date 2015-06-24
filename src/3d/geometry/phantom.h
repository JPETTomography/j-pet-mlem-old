#pragma once

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif


#include "3d/geometry/phantom_region.h"


namespace PET3D {


/*   Phantom -------------------------------------------------------------
 */

template <typename FType, typename SType, typename RNGType> class Phantom {
 public:
  using F = FType;
  using S = SType;
  using Event = PET3D::Event<F>;
  using RNG = RNGType;

 private:
  std::vector<PhantomRegion<F, RNG>*> region_list;
  std::vector<F> CDF;

  std::vector<Event> events;
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
#if DEBUG
      std::cerr << "CDF[" << i << "]=" << CDF[i];
#endif
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

  template <typename Generator> Event gen_event(Generator& generator) {
  again:
    size_t i_region = choose_region(generator);
    Point<F> p = region_list[i_region]->random_point(generator);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j]->in(p))
        goto again;
    }
    return Event(p, region_list[i_region]->random_direction(generator));
  }
};
}

#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>

#include "2d/geometry/event.h"
#include "2d/geometry/ellipse.h"
#include "2d/geometry/rectangle.h"
#include "2d/geometry/phantom_region.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif
const double RADIAN = M_PI / 180;
namespace PET2D {
/// Virtual phantom made of  regions
template <typename FType, typename SType> class Phantom {
 public:
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using RNG = std::minstd_rand0;
  using Event = PET2D::Event<F>;

 private:
  int n_events_;

  std::vector<PhantomRegion<F, RNG>*> region_list;
  std::vector<F> CDF;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom() : uniform_angle(-1, 1){};
  Phantom(const std::vector<PhantomRegion<F, RNG>*>& el)
      : uniform_angle(-1, 1), region_list(el) {

    calculate_cdf();
  }

  void calculate_cdf() {
    CDF.assign(region_list.size(), 0);
    CDF[0] = region_list[0]->weight();
    for (size_t i = 1; i < region_list.size(); i++) {
      CDF[i] = region_list[i]->weight() + CDF[i - 1];
    }
    F norm = CDF[region_list.size() - 1];
    for (size_t i = 0; i < region_list.size(); i++) {
      CDF[i] /= norm;
    }
  }

  void push_back_region(PhantomRegion<F, RNG>* region) {
    region_list.push_back(region);
  }

  int read_from_stream(std::istream& infile) {
    std::string line;
    int count = 0;
    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      std::string type;
      iss >> type;
      if (type == "ellipse") {
        double x, y, a, b, angle, acceptance;

        // on error
        if (!(iss >> x >> y >> a >> b >> angle >> acceptance))
          break;

        Ellipse<F> el(x, y, a, b, angle * RADIAN);

        auto region =
            new PET2D::EllipticalPhantomRegion<F, RNG>(el, acceptance);
        push_back_region(region);
        count++;
      } else if (type == "rectangle") {
        double x, y, a, b, acceptance;

        // on error
        if (!(iss >> x >> y >> a >> b >> acceptance))
          break;

        Rectangle<F> rec(x, y, a, b);

        auto region =
            new PET2D::RectangularPhantomRegion<F, RNG>(rec, acceptance);
        push_back_region(region);
        count++;
      } else {
        std::cerr << "unknow phantom type" << std::endl;
        exit(-1);
      }
    }
    return count;
  }

  template <class RNG> size_t choose_region(RNG& rng) {
    F r = uniform(rng);
    size_t i = 0;

    while (r > CDF[i])
      ++i;

    return i;
  }

  template <class RNG> Point<F> gen_point(RNG& rng) {
  again:
    size_t i_region = choose_region(rng);
    Point<F> p = region_list[i_region]->random_point(rng);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j]->contains(p))
        goto again;
    }
    return p;
  }

  template <class RNG> PET2D::Event<F> gen_event(RNG& rng) {
    Point<F> p = gen_point(rng);
    F rangle = F(M_PI_2) * uniform_angle(rng);
    return PET2D::Event<F>(p, rangle);
  }
};

}  // PET2D

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
const double RADIAN = M_PI / 180;
namespace PET2D {
namespace Strip {

// template <typename F> int sgn(F val) { return (0 < val) - (val < 0); }

/// Virtual phantom region made of ellipse and intensity
template <typename FType, typename RNG> struct PhantomRegion {
  using F = FType;
  using Point = PET2D::Point<F>;

  PhantomRegion(F intensity) : intensity(intensity) {}

  virtual bool contains(Point p) const = 0;
  virtual F weight() const = 0;
  virtual Point random_point(RNG&) = 0;

  const F intensity;
};

template <typename Shape, typename RNG>
class ShapePhantomRegion : public PhantomRegion<typename Shape::F, RNG> {
 public:
  using F = typename Shape::F;

  ShapePhantomRegion(const Shape& shape, F intensity)
      : PhantomRegion<F, RNG>(intensity),
        shape(shape),
        weight_(intensity * shape.area) {}

  bool contains(Point<F> p) const { return shape.contains(p); }

  const Shape shape;

  F weight() const { return weight_; }

 private:
  const F weight_;
};

template <typename FType, typename RNG>
class EllipticalPhantomRegion
    : public ShapePhantomRegion<PET2D::Ellipse<FType>, RNG> {
 public:
  using F = FType;
  using Ellipse = PET2D::Ellipse<F>;
  using Point = PET2D::Point<F>;

  EllipticalPhantomRegion(const Ellipse& ellipse, F intensity)
      : ShapePhantomRegion<Ellipse, RNG>(ellipse, intensity), gen_(ellipse) {}

  Point random_point(RNG& rng) { return gen_(rng); }

 private:
  EllipsePointGenerator<F> gen_;
};

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
            new PET2D::Strip::EllipticalPhantomRegion<F, RNG>(el, acceptance);
        push_back_region(region);
      } else {
        std::cerr << "unknow phantom type" << std::endl;
        exit(-1);
      }
    }
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
    Point<F> p = region_list[i_region]->random_point(generator);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j]->contains(p))
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

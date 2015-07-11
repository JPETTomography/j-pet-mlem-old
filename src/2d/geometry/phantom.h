#pragma once

#include <cmath>
#include <vector>
#include <ctime>
#include <algorithm>

#include "2d/geometry/event.h"
#include "2d/geometry/ellipse.h"
#include "2d/geometry/rectangle.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

const double RADIAN = M_PI / 180;

namespace PET2D {

/// Virtual phantom made of regions
template <class RNGClass, typename FType> class Phantom {
 public:
  using RNG = RNGClass;
  using F = FType;
  using Point = PET2D::Point<F>;
  using Event = PET2D::Event<F>;

  /// Abstract phantom region (must subclass)

  /// Must provide at least intensity for the region.
  struct Region {
    Region(F intensity) : intensity(intensity) {}

    virtual bool contains(Point p) const = 0;
    virtual F weight() const = 0;
    virtual Point random_point(RNG&) = 0;

    const F intensity;
  };

  /// Region that represent custom shape
  template <class Shape> class ShapeRegion : public Region {
   public:
    ShapeRegion(const Shape& shape, F intensity)
        : Region(intensity), shape(shape), weight_(intensity * shape.area) {}

    bool contains(Point p) const { return shape.contains(p); }

    const Shape shape;

    F weight() const { return weight_; }

   private:
    const F weight_;
  };

  /// Elliptical region
  class EllipticalRegion : public ShapeRegion<PET2D::Ellipse<F>> {
   public:
    using Ellipse = PET2D::Ellipse<F>;

    EllipticalRegion(const Ellipse& ellipse, F intensity)
        : ShapeRegion<Ellipse>(ellipse, intensity), gen_(ellipse) {}

    Point random_point(RNG& rng) { return gen_(rng); }

   private:
    EllipsePointGenerator<F> gen_;
  };

  /// Rectangular region
  class RectangularRegion : public ShapeRegion<PET2D::Rectangle<F>> {
   public:
    using Rectangle = PET2D::Rectangle<F>;

    RectangularRegion(const Rectangle& rectangle, F intensity)
        : ShapeRegion<Rectangle>(rectangle, intensity), gen_(rectangle) {}

    Point random_point(RNG& rng) { return gen_(rng); }

   private:
    PET2D::RectanglePointGenerator<F> gen_;
  };

  using RegionPtrList = std::vector<Region*>;

 private:
  int n_events_;

  RegionPtrList region_list;
  std::vector<F> CDF;

  std::uniform_real_distribution<F> uniform;
  std::uniform_real_distribution<F> uniform_angle;

 public:
  Phantom() : uniform_angle(-1, 1) {}
  Phantom(const RegionPtrList& el) : uniform_angle(-1, 1), region_list(el) {

    calculate_cdf();
  }

  void calculate_cdf() {
    if (region_list.size() == 0) {
      throw("must specify at least one region");
    }
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

  void push_back_region(Region* region) { region_list.push_back(region); }

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

        auto region = new EllipticalRegion(el, acceptance);
        push_back_region(region);
        count++;
      } else if (type == "rectangle") {
        double x, y, a, b, acceptance;

        // on error
        if (!(iss >> x >> y >> a >> b >> acceptance))
          break;

        Rectangle<F> rec(x, y, a, b);

        auto region = new RectangularRegion(rec, acceptance);
        push_back_region(region);
        count++;
      } else {
        std::cerr << "unknow phantom type" << std::endl;
        exit(-1);
      }
    }
    return count;
  }

  size_t choose_region(RNG& rng) {
    F r = uniform(rng);
    size_t i = 0;

    while (r > CDF[i])
      ++i;

    return i;
  }

  Point gen_point(RNG& rng) {
  again:
    size_t i_region = choose_region(rng);
    Point p = region_list[i_region]->random_point(rng);
    for (size_t j = 0; j < i_region; j++) {
      if (region_list[j]->contains(p))
        goto again;
    }
    return p;
  }

  Event gen_event(RNG& rng) {
    Point p = gen_point(rng);
    F rangle = F(M_PI_2) * uniform_angle(rng);
    return Event(p, rangle);
  }
};

}  // PET2D

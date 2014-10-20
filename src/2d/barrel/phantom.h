#pragma once

#include <cmath>
#include <vector>

#include "2d/geometry/point.h"

namespace PET2D {
namespace Barrel {

/// Elliptical phantom region
template <typename FType = double> struct EllipticalRegion {
  typedef FType F;
  typedef PET2D::Point<F> Point;

  const Point center;
  const F a;
  const F b;
  const F phi;
  const F activity;

  EllipticalRegion(Point center, F a, F b, F phi, F activity)
      : center(center), a(a), b(b), phi(phi), activity(activity) {
    sincos(phi, &sin, &cos);
    inv_a2 = 1 / (a * a);
    inv_b2 = 1 / (b * b);
  }

#if !__CUDACC__
  EllipticalRegion(std::istream& in)
      : EllipticalRegion(Point(in),
                         read<F>(in),
                         read<F>(in),
                         read<F>(in),
                         read<F>(in)) {}
#endif

  bool operator()(Point p) const {

    auto d = p - center;

    auto r_x = d.x * cos + d.y * sin;
    auto r_y = -d.x * sin + d.y * cos;

    auto r2 = r_x * r_x * inv_a2 + r_y * r_y * inv_b2;

    return r2 < static_cast<F>(1);
  }

 private:
  F inv_a2;
  F inv_b2;
  F sin;
  F cos;

#if __APPLE__ || _MSC_VER
  template <typename F> static inline void sincos(F a, F* s, F* c) {
    *s = std::sin(a);
    *c = std::cos(a);
  }
#endif
};

/// Phantom made of elliptical regions
template <typename FType = double>
class Phantom : public std::vector<EllipticalRegion<FType>> {
 public:
  typedef FType F;
  typedef PET2D::Point<F> Point;

  size_t n_regions() const { return this->size(); }

  F activity(Point p) const {
    for (auto& region : *this) {
      if (region(p))
        return region.activity;
    }
    return static_cast<F>(0);
  }

  bool test_emit(Point p, F rnd) const { return activity(p) > rnd; }
};

/// Single point source
template <typename FType = double> struct PointSource {
  typedef FType F;
  typedef PET2D::Point<F> Point;

  PointSource(Point p, F intensity) : p(p), intensity(intensity) {}

  PointSource(std::istream& in)
      : p(read<F>(in), read<F>(in)), intensity(read<F>(in)) {}

  const Point p;
  const F intensity;
};

/// Phantom made of point sources
template <typename FType = double>
class PointSources : public std::vector<PointSource<FType>> {
 public:
  typedef FType F;

  size_t n_sources() const { return this->size(); }

  void normalize() {
    total = static_cast<F>(0);
    for (auto& source : *this) {
      total += source.intensity;
    }

    F cumulant = static_cast<F>(0);
    cumulants.clear();
    for (auto& source : *this) {
      cumulant += source.intensity / total;
      cumulants.push_back(cumulant);
    }
  }

  Point<F> draw(F rng) {
    int i = 0;
    while (rng > cumulants[i]) {
      ++i;
    }
    return (*this)[i].p;
  }

 private:
  std::vector<F> cumulants;
  F total;
};
}  // Barrel
}  // PET2D

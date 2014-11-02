#pragma once

#include <cmath>
#include <vector>

#include "2d/geometry/ellipse.h"

namespace PET2D {
namespace Barrel {

/// Virtual phantom made of elliptical regions
template <typename FType = double>
class Phantom : public std::vector<EllipticalSource<FType>> {
 public:
  typedef FType F;
  typedef PET2D::Point<F> Point;

  size_t n_regions() const { return this->size(); }

  F intensity(Point p) const {
    for (auto& region : *this) {
      if (region.contains(p))
        return region.intensity;
    }
    return 0;
  }

  bool test_emit(Point p, F rnd) const { return intensity(p) > rnd; }
};

/// Virtual phantom made of point sources
template <typename FType = double>
class PointPhantom : public std::vector<PointSource<FType>> {
 public:
  typedef FType F;

  size_t n_sources() const { return this->size(); }

  void normalize() {
    total = 0;
    for (auto& source : *this) {
      total += source.intensity;
    }

    F cumulant = 0;
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
    return (*this)[i];
  }

 private:
  std::vector<F> cumulants;
  F total;
};
}  // Barrel
}  // PET2D

#pragma once

#include <cmath>
#include <vector>
#include "2d/geometry/point.h"

template <typename FType = double> struct EllipticalRegion {
  typedef FType F;
  typedef ::Point<F> Point;

  EllipticalRegion(Point center, F a, F b, F phi, F activity)
      : center_(center), a_(a), b_(b), phi_(phi), activity_(activity) {
    init();
  }

  EllipticalRegion(std::istream& in) {
    in >> center_.x >> center_.y >> a_ >> b_ >> phi_ >> activity_;
    init();
  }

  F activity() const { return activity_; }

  bool operator()(Point p) const {

    auto d = p - center_;

    auto r_x = d.x * cos_ + d.y * sin_;
    auto r_y = -d.x * sin_ + d.y * cos_;

    auto r2 = r_x * r_x * inv_a2 + r_y * r_y * inv_b2;

    return r2 < static_cast<F>(1);
  }

  F center() const { return center_; }
  F a() const { return a_; }
  F b() const { return b_; }
  F phi() const { return phi_; }

 private:
  void init() {
    sincos(phi_, &sin_, &cos_);
    inv_a2 = static_cast<F>(1) / (a_ * a_);
    inv_b2 = static_cast<F>(1) / (b_ * b_);
  }

  Point center_;
  F a_;
  F b_;

  F phi_;
  F activity_;

  F inv_a2;
  F inv_b2;

  F sin_;
  F cos_;

#if __APPLE__ || _MSC_VER
  template <typename F> static inline void sincos(F a, F* s, F* c) {
    *s = std::sin(a);
    *c = std::cos(a);
  }
#endif
};

template <typename FType = double>
class Phantom : public std::vector<EllipticalRegion<FType>> {
 public:
  typedef FType F;
  typedef ::Point<F> Point;

  size_t n_regions() const { return this->size(); }

  F activity(Point p) const {
    for (auto& region : *this) {
      if (region(p))
        return region.activity();
    }
    return static_cast<F>(0);
  }

  bool test_emit(Point p, F rnd) const { return activity(p) > rnd; }
};

template <typename FType = double> struct PointSource {
  typedef FType F;
  typedef ::Point<F> Point;

  PointSource(Point p, F intensity) : p(p), intensity(intensity) {}

  PointSource(std::istream& in) { in >> p.x >> p.y >> intensity; }

  Point p;
  F intensity;
};

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
      cumulant += source.intensity /= total;
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

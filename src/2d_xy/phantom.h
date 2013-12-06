#pragma once

#include <cmath>
#include <vector>
#include "geometry/point.h"

class EllipticalRegion {
 public:
  EllipticalRegion(double x,
                   double y,
                   double a,
                   double b,
                   double phi,
                   double act)
      : x_(x), y_(y), a_(a), b_(b), phi_(phi), activity_(act) {
    sincos(phi, &sin_, &cos_);
    inv_a2 = 1.0 / (a_ * a_);
    inv_b2 = 1.0 / (b_ * b_);
  }

  double activity() const { return activity_; }

  bool in(double x, double y) const {

    double dx = x - x_;
    double dy = y - y_;

    double r_x = dx * cos_ + dy * sin_;
    double r_y = -dx * sin_ + dy * cos_;

    double r2 = r_x * r_x * inv_a2 + r_y * r_y * inv_b2;

    return r2 < 1.0;
  }

  double x() const { return x_; }
  double y() const { return y_; }
  double a() const { return a_; }
  double b() const { return b_; }
  double phi() const { return phi_; }

 private:
  double x_;
  double y_;
  double a_;
  double b_;

  double phi_;
  double activity_;

  double inv_a2;
  double inv_b2;

  double sin_;
  double cos_;

#ifdef __APPLE__
  template <typename T> static inline void sincos(T a, T* s, T* c) {
    *s = std::sin(a);
    *c = std::cos(a);
  }
#endif
};

class Phantom : public std::vector<EllipticalRegion> {

 public:
  size_t n_regions() const { return size(); }

  void add_region(double x,
                  double y,
                  double a,
                  double b,
                  double phi,
                  double act) {
    this->push_back(EllipticalRegion(x, y, a, b, phi, act));
  }

  double activity(double x, double y) const {
    for (auto rit = this->rbegin(); rit != this->rend(); ++rit) {
      if (rit->in(x, y)) {
        return rit->activity();
      }
    }
    return 0.0;
  }

  bool test_emit(double x, double y, double rnd) const {
    return activity(x, y) > rnd;
  }
};

template <typename FType = double> struct PointSource {
  typedef FType F;

  PointSource(F x, F y, F intensity_a) : p(x, y), intensity(intensity_a) {}
  Point<F> p;
  F intensity;
};

template <typename FType = double>
class PointSources : public std::vector<PointSource<FType>> {
 public:
  typedef FType F;

  size_t n_sources() const { return this->size(); }

  void add(F x, F y, F intensity) {
    this->push_back(PointSource<F>(x, y, intensity));
  }

  void normalize() {
    total = 0.0;
    for (auto& source : *this) {
      total += source.intensity;
    }

    F cumulant = 0.0;
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

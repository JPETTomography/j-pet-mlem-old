#pragma once

#include "polygon.h"

template <typename F = double> class Detector : public Polygon<F> {
 public:
  typedef F Angle;
  typedef typename Polygon<F>::Point Point;

  Detector(F w, F h) {
    this->push_back(Point(w / 2., h / 2.));
    this->push_back(Point(w / 2., -h / 2.));
    this->push_back(Point(-w / 2., -h / 2.));
    this->push_back(Point(-w / 2., h / 2.));
  }

  Detector rotated(Angle phi) {
    Detector r;
    for (auto it = this->begin(); it != this->end(); ++it) {
      r.push_back((*it).rotated(phi));
    }
    return r;
  }

  Detector& operator+=(Point t) {
    for (auto it = this->begin(); it != this->end(); ++it)
      *it += t;
    return *this;
  }

 private:
  Detector() {}
};

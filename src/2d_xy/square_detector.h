#pragma once

#include "geometry/polygon.h"

template <typename FType = double>
class SquareDetector : public Polygon<4, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<4, F>::Point Point;

  SquareDetector(F w, F h) {
    this->push_back(Point(w / 2., h / 2.));
    this->push_back(Point(w / 2., -h / 2.));
    this->push_back(Point(-w / 2., -h / 2.));
    this->push_back(Point(-w / 2., h / 2.));
  }

  SquareDetector& rotate(Angle phi) {
    SquareDetector r;
    for (auto it = this->begin(); it != this->end(); ++it) {
      it->rotate(phi);
    }
    return *this;
  }

  SquareDetector rotated(Angle phi) {
    SquareDetector r;
    for (auto it = this->begin(); it != this->end(); ++it) {
      r.push_back((*it).rotated(phi));
    }
    return r;
  }

  SquareDetector& operator+=(Point t) {
    for (auto it = this->begin(); it != this->end(); ++it)
      *it += t;
    return *this;
  }

 private:
  SquareDetector() {}
};

#pragma once

#include "geometry/polygon.h"

template <typename FType = double>
class TriangleDetector : public Polygon<4, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<3, F>::Point Point;

  TriangleDetector(F w, F h) {
    this->push_back(Point(w / 2., h / 2.));
    this->push_back(Point(-w / 2., -h / 2.));
    this->push_back(Point(w / 2., h));
  }

  TriangleDetector& rotate(Angle phi) {
    TriangleDetector r;
    for (auto& p : *this) {
      p.rotate(phi);
    }
    return *this;
  }

  TriangleDetector rotated(Angle phi) {
    TriangleDetector r;
    for (auto& p : *this) {
      r.push_back(p.rotated(phi));
    }
    return r;
  }

  TriangleDetector& operator+=(Point t) {
    for (auto& p : *this) {
      p += t;
    }
    return *this;
  }

 private:
  TriangleDetector() {}
};

#pragma once

#include <cmath>
#include "geometry/polygon.h"

template <typename FType = double>
class HexagonalTriangle : public Polygon<6, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<6, F>::Point Point;

  HexagonalTriangle(F w, F h, F size) {

    for(int i = 0; i < 6; ++i){

      angle = (2*M_PI)/F(6.0) * (i + F(0.5));

      this->push_back(Point(w + size * std::cos(angle)  , h + size * std::sin(angle) ));

    }
  }

  HexagonalTriangle& rotate(Angle phi) {
    TriangleDetector r;
    for (auto& p : *this) {
      p.rotate(phi);
    }
    return *this;
  }

  HexagonalTriangle rotated(Angle phi) {
    TriangleDetector r;
    for (auto& p : *this) {
      r.push_back(p.rotated(phi));
    }
    return r;
  }

  HexagonalTriangle& operator+=(Point t) {
    for (auto& p : *this) {
      p += t;
    }
    return *this;
  }

 private:
  HexagonalTriangle() {}
};

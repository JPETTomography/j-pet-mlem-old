#pragma once

#include "geometry/polygon.h"

template <typename FType = double>
class SquareDetector : public Polygon<4, FType> {
 public:
  typedef FType F;
  typedef typename Polygon<4, F>::Point Point;

  SquareDetector(F w, F h, F d __attribute__((unused)) = F()) {
    this->push_back(Point(w / 2, h / 2));
    this->push_back(Point(w / 2, -h / 2));
    this->push_back(Point(-w / 2, -h / 2));
    this->push_back(Point(-w / 2, h / 2));
  }

  static F default_height_for_width(const F w) { return w; }

 private:
  SquareDetector() {}
};

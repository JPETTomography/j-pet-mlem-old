#pragma once

#include "2d/geometry/polygon.h"

namespace PET2D {
namespace Barrel {

template <typename FType = double>
class SquareDetector : public Polygon<4, FType> {
 public:
  typedef FType F;
  typedef typename Polygon<4, F>::Point Point;

  SquareDetector(F w, F h, F d) {
    (void)d;  // unused
    this->push_back(Point(w / 2, h / 2));
    this->push_back(Point(w / 2, -h / 2));
    this->push_back(Point(-w / 2, -h / 2));
    this->push_back(Point(-w / 2, h / 2));
  }

  static F default_height_for_width(const F w) { return w; }

 private:
  SquareDetector() {}
};
}  // Barrel
}  // PET2D

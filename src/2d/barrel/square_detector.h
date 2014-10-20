#pragma once

#include "2d/geometry/polygon.h"

namespace PET2D {
namespace Barrel {

/// Single square detector
template <typename FType = double>
class SquareDetector : public Polygon<4, FType> {
 public:
  typedef FType F;
  typedef typename Polygon<4, F>::Point Point;

  SquareDetector(F w, F h, F d) {
    (void)d;  // unused
    this->emplace_back(w / 2, h / 2);
    this->emplace_back(w / 2, -h / 2);
    this->emplace_back(-w / 2, -h / 2);
    this->emplace_back(-w / 2, h / 2);
  }

  static F default_height_for_width(const F w) { return w; }

 private:
  SquareDetector() {}
};
}  // Barrel
}  // PET2D

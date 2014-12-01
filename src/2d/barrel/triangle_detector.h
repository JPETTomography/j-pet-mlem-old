#pragma once

#include "polygonal_detector.h"

namespace PET2D {
namespace Barrel {

/// Single triangular detector
template <typename FType = double>
class TriangleDetector : public PolygonalDetector<3, FType> {
 public:
  using F = FType;
  using Angle = F;
  using Point = typename Polygon<3, F>::Point;

  TriangleDetector(F w, F h, F d = F()) {
    if (d > F()) {
      this->emplace_back(w / 2, d / 2 - h);
      this->emplace_back(-w / 2, d / 2 - h);
      this->emplace_back(0, d / 2);
    } else {
      this->emplace_back(w / 2, -h / 2);
      this->emplace_back(-w / 2, -h / 2);
      this->emplace_back(0, h / 2);
    }
  }

  static F default_height_for_width(const F w) {
    return w * std::sqrt(static_cast<F>(3)) / 2;
  }

 private:
  TriangleDetector() {}
};
}  // Barrel
}  // PET2D

#pragma once

#include "polygonal_detector.h"

namespace PET2D {
namespace Barrel {

/// Single triangular detector
////
/// Represents single detector with triangular shape:
/// \image html shape_triangle.pdf.png
template <typename FType>
class TriangleDetector : public PolygonalDetector<3, FType> {
 public:
  using Base = PolygonalDetector<3, FType>;
  using F = FType;
  using Angle = F;
  using Point = typename Polygon<3, F>::Point;

  TriangleDetector(F w, F h, F d = 0) {
    if (d > 0) {
      this->emplace_back(w / 2, d / 2 - h);
      this->emplace_back(-w / 2, d / 2 - h);
      this->emplace_back(0, d / 2);
    } else {
      this->emplace_back(w / 2, -h / 2);
      this->emplace_back(-w / 2, -h / 2);
      this->emplace_back(0, h / 2);
    }
  }

  TriangleDetector(Base&& base) : Base(std::forward<Base>(base)) {}
#if !_MSC_VER
  TriangleDetector(typename Base::Base&& base)
      : Base(std::forward<typename Base::Base>(base)) {}
#endif

  static F default_height_for_width(const F w) {
    return w * std::sqrt(static_cast<F>(3)) / 2;
  }

 private:
  TriangleDetector() {}
};
}  // Barrel
}  // PET2D

#pragma once

#include "polygonal_detector.h"

namespace PET2D {
namespace Barrel {

/// Single square detector
////
/// Represents single detector with rectangular shape:
/// \image html shape_rectangle.pdf.png
template <typename FType>
class SquareDetector : public PolygonalDetector<4, FType> {
 public:
  using Base = PolygonalDetector<4, FType>;
  using F = FType;
  using Point = typename Polygon<4, F>::Point;

  SquareDetector(F w, F h, F d) {
    (void)d;  // unused
    this->emplace_back(w / 2, h / 2);
    this->emplace_back(w / 2, -h / 2);
    this->emplace_back(-w / 2, -h / 2);
    this->emplace_back(-w / 2, h / 2);
  }

  SquareDetector(F x, F y, F w, F h) {
    this->emplace_back(x + w / 2, y + h / 2);
    this->emplace_back(x + w / 2, y - h / 2);
    this->emplace_back(x - w / 2, y - h / 2);
    this->emplace_back(x - w / 2, y + h / 2);
  }

  SquareDetector(Base&& base) : Base(std::forward<Base>(base)) {}
  SquareDetector(const Base& base) : Base(base) {}
#if !_MSC_VER
  SquareDetector(typename Base::Base&& base)
      : Base(std::forward<typename Base::Base>(base)) {}
#endif

  F width() const { return ((*this)[0] - (*this)[3]).length(); }
  F height() const { return ((*this)[0] - (*this)[1]).length(); }

  static F default_height_for_width(const F w) { return w; }

 private:
  SquareDetector() {}
};
}  // Barrel
}  // PET2D

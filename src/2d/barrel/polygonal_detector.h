#pragma once

#include "2d/geometry/polygon.h"

namespace PET2D {
namespace Barrel {

/// Single detector with shape of custom polygon
template <std::size_t NVertices, typename FType = double>
class PolygonalDetector : public Polygon<NVertices, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<NVertices, F>::Point Point;

  PolygonalDetector(F w, F h, F d = F()) {
    (void)h, (void)d;  // unused

    auto radius = w / (4 * std::sin(F(M_PI) / F(NVertices)));
    auto step = 2 * F(M_PI) / F(NVertices);

    for (std::size_t i = 0; i < NVertices; ++i) {
      auto angle = step * (i + F(0.5));
      this->emplace_back(radius * std::cos(angle), radius * std::sin(angle));
    }
  }

  static F default_height_for_width(const F w) { return w; }

 private:
  PolygonalDetector() {}
};
}  // Barrel
}  // PET2D

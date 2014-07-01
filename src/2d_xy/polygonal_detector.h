#pragma once

#include "geometry/polygon.h"

template <std::size_t NVertices, typename FType = double>
class PolygonalDetector : public Polygon<NVertices, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<NVertices, F>::Point Point;

  PolygonalDetector(F w,
                    F h __attribute__((unused)),
                    F d __attribute__((unused)) = F()) {
    auto radius = w / (4 * std::sin(F(M_PI) / F(NVertices)));
    auto step = 2 * F(M_PI) / F(NVertices);

    for (std::size_t i = 0; i < NVertices; ++i) {
      auto angle = step * (i + F(0.5));
      this->push_back(
          Point(radius * std::cos(angle), radius * std::sin(angle)));
    }
  }

  static F default_height_for_width(const F w) { return w; }

 private:
  PolygonalDetector() {}
};

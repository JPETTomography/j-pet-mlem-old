#pragma once

#include "geometry/polygon.h"

template <std::size_t NVertices, typename FType = double>
class PolygonalDetector : public Polygon<NVertices, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<NVertices, F>::Point Point;

  PolygonalDetector(F w, F h __attribute__((unused))) {
    auto radius = w / 2.;
    auto step = 2. * M_PI / F(NVertices);

    for (std::size_t i = 0; i < NVertices; ++i) {
      auto angle = step * (i + F(0.5));
      this->push_back(
          Point(radius * std::cos(angle), radius * std::sin(angle)));
    }
  }

  static F default_height_for_width(const F w __attribute__((unused))) {
    return F();
  }

 private:
  PolygonalDetector() {}
};

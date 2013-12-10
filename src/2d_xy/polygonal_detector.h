#pragma once

#include "geometry/polygon.h"

template <std::size_t NVertices, typename FType = double>
class PolygonalDetector : public Polygon<NVertices, FType> {
 public:
  typedef FType F;
  typedef F Angle;
  typedef typename Polygon<NVertices, F>::Point Point;

  PolygonalDetector(F w, F h, F size) {
    for (auto i = 0; i < NVertices; ++i) {
      angle = (2. * M_PI) / F(NVertices) * (i + F(0.5));

      this->push_back(
          Point(w + size * std::cos(angle), h + size * std::sin(angle)));
    }
  }

 private:
  PolygonalDetector() {}
};

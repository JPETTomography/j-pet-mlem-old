#pragma once

#include "circle_detector.h"
#include "2d/geometry/polygon.h"

#if !__CUDACC__
#include "util/json.h"
#endif

namespace PET2D {
namespace Barrel {

/// Single detector with shape of custom polygon
////
/// Represents single detector with convex polygonal shape, such as hexagon:
/// \image html shape_hexagon.pdf.png
template <std::size_t NVertices, typename FType>
class PolygonalDetector : public Polygon<NVertices, FType> {
 public:
  static std::size_t const N = NVertices;
  using Base = Polygon<NVertices, FType>;
  using F = FType;
  using Angle = F;
  using Point = typename Polygon<NVertices, F>::Point;
  using CircleDetector = Barrel::CircleDetector<F>;

  PolygonalDetector(F w, F h, F d = 0) {
    (void)h, (void)d;  // unused

    auto radius = w / (4 * compat::sin(F(M_PI) / F(NVertices)));
    auto step = 2 * F(M_PI) / F(NVertices);

    for (std::size_t i = 0; i < NVertices; ++i) {
      auto angle = step * (i + F(0.5));
      this->emplace_back(radius * compat::cos(angle),
                         radius * compat::sin(angle));
    }
  }

  PolygonalDetector(Base&& base) : Base(std::forward<Base>(base)) {}
  PolygonalDetector(const Base& rhs) : Base(rhs) {}

  /// \returns circumscribed circular detector
  CircleDetector circumscribe_circle() const {
    Point center = this->center();
    F radius = 0;
    for (const Point& p : *this) {
      radius = compat::max(radius, (p - center).length());
    }
    return CircleDetector(radius, center);
  }

  static F default_height_for_width(const F w) { return w; }

 protected:
  PolygonalDetector() {}
};
}  // Barrel
}  // PET2D

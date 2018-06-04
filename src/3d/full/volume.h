#ifndef VOLUME_H
#define VOLUME_H

#include "3d/geometry/ray.h"

namespace PET3D {
namespace Full {

template <typename F> class Material {};

template <typename F> class Volume {
 public:
  virtual ray_tracing::intersection_result<F> intersects_with(
      const ray_tracing::Ray<F>& ray) = 0;
  bool is_detector;

 private:
  Material<F> material;
};

template <typename F> class BoxVolume : public Volume<F> {
 public:
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;
  BoxVolume(const Point& center,
            const Vector& a_u,
            const Vector& a_v,
            const Vector& a_w,
            F h_u,
            F h_v,
            F h_w)
      : box(center, a_u, a_v, a_w, h_u, h_v, h_w) {}

  virtual ray_tracing::intersection_result<F> intersects_with(
      const ray_tracing::Ray<F>& ray) {
    return ray_tracing::intersect(ray, box);
  }

  static BoxVolume* AAB(const Point& p1, const Point& p2) {
    Vector half_diag = F(0.5) * (p2 - p1);
    return new BoxVolume(p1 + half_diag,
                         Vector::e_x(),
                         Vector::e_y(),
                         Vector::e_z(),
                         half_diag.x,
                         half_diag.y,
                         half_diag.z);
  }

 private:
  ray_tracing::Box<F> box;
};
}
}
#endif  // VOLUME_H

#ifndef RAY_H
#define RAY_H

#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"

namespace ray_tracing {
template <typename F> class Ray {
 public:
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;

  Ray(const Point& p, const Vector& d) : p(p), d(d) {}
  const Point p;
  const Vector d;
};

template <typename F> class Box {
 public:
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;

  Box(const Point& center,
      const Vector& a_u,
      const Vector& a_v,
      const Vector& a_w,
      F h_u,
      F h_v,
      F h_w)
      : center(center),
        a_u(a_u),
        a_v(a_v),
        a_w(a_w),
        h_u(h_u),
        h_v(h_v),
        h_w(h_w) {}

 public:
  const Point center;
  const Vector a_u, a_v, a_w;
  const F h_u, h_v, h_w;

  static Box AAB(const Point& p1, const Point& p2) {
    Vector half_diag = F(0.5) * (p2 - p1);
    return Box(p1 + half_diag,
               Vector::e_x(),
               Vector::e_y(),
               Vector::e_z(),
               half_diag.x,
               half_diag.y,
               half_diag.z);
  }
};

template <typename F> using intersection_result = std::pair<bool, F>;

template <typename F>
intersection_result<F> intersect(const Ray<F>& ray, const Box<F>& box) {
  return std::make_pair(true, 0.5);
}
}

#endif  // RAY_H

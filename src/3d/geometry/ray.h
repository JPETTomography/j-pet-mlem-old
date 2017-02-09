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

  Point operator()(F t) { return p + t * d; }
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
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;

  F t_min = std::numeric_limits<F>::lowest();
  F t_max = std::numeric_limits<F>::max();

  Vector p = box.center - ray.p;

  Vector normals[3] = { box.a_u, box.a_v, box.a_w };
  F sides[3] = { box.h_u, box.h_v, box.h_w };

  for (int i = 0; i < 3; i++) {
    Vector a = normals[i];
    F e = a.dot(p);
    F f = a.dot(ray.d);
    if (std::fabs(f) > std::numeric_limits<F>::epsilon()) {
      F t1 = (e + sides[i]) / f;
      F t2 = (e - sides[i]) / f;
      if (t1 > t2)
        std::swap(t1, t2);
      if (t1 > t_min)
        t_min = t1;
      if (t2 < t_max)
        t_max = t2;
      if (t_min > t_max)
        return std::make_pair(false, 0);
      if (t_max < 0)
        return std::make_pair(false, 0);
    } else {

      if ((-e - sides[i]) > 0 || (-e + sides[i]) < 0)
        return std::make_pair(false, 0);
    }
  }

  if (t_min > 0)
    return std::make_pair(true, t_min);
  else
    return std::make_pair(true, t_max);
}
}

#endif  // RAY_H

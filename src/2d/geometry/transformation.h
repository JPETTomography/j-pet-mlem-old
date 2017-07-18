#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "vector.h"
#include "point.h"
namespace PET2D {
template <typename FType> class Transformation {
 public:
  using F = FType;
  using Vector = PET2D::Vector<F>;
  using Point = PET2D::Point<F>;

  Transformation(F r, Vector t, bool mirror_x = false)
      : rotation(r), translation(t), mirror_x(mirror_x){};
  Transformation(F r) : Transformation(r, Vector(0, 0)){};
  Transformation(Vector t) : Transformation(0, t){};
  Transformation() : Transformation(0, Vector(0, 0)){};

  Point operator()(Point p) {
    if (mirror_x) {
      auto q = Point(-p.x, p.y);
      return q.rotated(rotation) + translation;
    } else {
      p.rotated(rotation) + translation;
    }
  }
  Vector operator()(Vector v) { return v.rotated(rotation); }

  const F rotation;
  const Vector translation;
  const bool mirror_x;
};

template <typename F>
Transformation<F> operator*(const Transformation<F>& t2,
                            const Transformation<F>& t1) {
  return Transformation<F>(
      t2.rotation + t1.rotation,
      t1.translation.rotated(t2.rotation) + t2.translation);
}
}
#endif  // TRANSFORMATION_H

#ifndef KERNEL_TESTING_TOOLBOX
#define KERNEL_TESTING_TOOLBOX

#include <cmath>

#include "3d/geometry/vector.h"

namespace PET2D {

namespace Strip {
namespace Testing {

template <typename F> using Vector3D = PET3D::Vector<F>;

template <typename F> class FrameEvent;

template <typename F> struct Event {

  Event(F x, F y, F theta) : tan(std::tan(theta)), theta(theta), y(y), x(x){};
  Event(const FrameEvent<F>& fe, F R)
      : tan((fe.zup - fe.zdn) / (F(2.0) * R)),
        theta(atan(tan)),
        y(-F(0.5) * fe.dl * std::cos(theta)),
        x(F(0.5) * (fe.zup + fe.zdn + 2 * y * tan)){};

  const F tan;
  const F theta;
  const F y;
  const F x;
};

template <typename F> struct FrameEvent {

  FrameEvent(F zup, F zdn, F dl) : zup(zup), zdn(zdn), dl(dl) {}
  FrameEvent(const Event<F> evt, F R)
      : zup(evt.x + (R - evt.y) * evt.tan),
        zdn(evt.x - (R + evt.y) * evt.tan),
        dl(-F(2) * evt.y * std::sqrt(F(1) + evt.tan * evt.tan)){};

  const F zup, zdn, dl;
};

template <typename F>
Vector3D<F> operator-(const FrameEvent<F>& fel, const FrameEvent<F>& fer) {
  return Vector3D<F>(fel.zup - fer.zup, fel.zdn - fer.zdn, fel.dl - fer.dl);
}

template <typename F>
F diagonal_product(const Vector3D<F>& diag, const Vector3D<F>& vec) {
  Vector3D<F> res(diag);
  res *= vec;
  return res.dot(vec);
}

template <typename F> F gauss(const Vector3D<F>& diag, const Vector3D<F>& vec) {
  return std::exp(-F(0.5) * diagonal_product(diag, vec));
}

template <typename F> F sensitivity(const FrameEvent<F>& fe, F L) {
  F l2 = L / 2;

  return (fe.zup <= l2 && fe.zup >= -l2 && fe.zdn <= l2 && fe.zdn >= -l2) ? 1.0
                                                                          : 0.0;
}
}
}
}

#endif  // KERNEL_TESTING_TOOLBOX
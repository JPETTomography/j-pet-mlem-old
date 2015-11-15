#ifndef KERNEL_TESTING_TOOLBOX
#define KERNEL_TESTING_TOOLBOX

#include <cmath>

#include "3d/geometry/vector.h"

namespace PET2D {

namespace Strip {
namespace Testing {

const double PI_FACTOR = 1 / std::pow(2 * M_PI, 3.0 / 2.0);

template <typename F> using Vector3D = PET3D::Vector<F>;

template <typename F> struct FrameEvent;

template <typename F> struct Event {

  Event(F x, F y, F theta) : tan(std::tan(theta)), theta(theta), y(y), x(x) {}
  Event(const FrameEvent<F>& fe, F R)
      : tan((fe.zup - fe.zdn) / (F(2.0) * R)),
        theta(atan(tan)),
        y(-F(0.5) * fe.dl * std::cos(theta)),
        x(F(0.5) * (fe.zup + fe.zdn + 2 * y * tan)) {}

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
        dl(-F(2) * evt.y * std::sqrt(F(1) + evt.tan * evt.tan)) {}

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

template <typename F>
F weight(const Vector3D<F>& diag,
         const FrameEvent<F>& meas,
         const FrameEvent<F>& exact,
         F L) {

  if (sensitivity(exact, L) > 0) {
    auto diff = meas - exact;
    auto g = gauss(diag, diff);
    auto det = diag.x * diag.y * diag.z;

    return PI_FACTOR * std::sqrt(det) * g * sensitivity(exact, L);
  }
  return 0;
}

template <typename F> F sensitivity(const FrameEvent<F>& fe, F L) {
  F l2 = L / 2;

  return (fe.zup <= l2 && fe.zup >= -l2 && fe.zdn <= l2 && fe.zdn >= -l2) ? 1.0
                                                                          : 0.0;
}

template <typename F> std::pair<F, F> theta_min_max(F x, F y, F R, F L) {
  const F l2 = 0.5 * L;
  return std::make_pair(
      std::atan(std::max(-(l2 + x) / (R - y), (-l2 + x) / (R + y))),
      std::atan(std::min((l2 - x) / (R - y), (l2 + x) / (R + y))));
}

template <typename F> F sensitivity(F x, F y, F R, F L) {
  auto theta = theta_min_max(x, y, R, L);

  return (theta.second - theta.first) / M_PI;
}

F theta_integral(const Vector3D<F> diag,
                 const FrameEvent<F>& evt,
                 F x,
                 F y,
                 F R,
                 F L,
                 F d = 0.01) {
  double sum = 0.0;

  auto theta_lim = theta_min_max(x, y, R, L);

  for (F theta = theta_lim.first; theta < theta_lim.second; theta += d) {
    FrameEvent<F> exact(Event<F>(x, y, theta), R);
    sum += weight(diag, evt, exact, L);
  }
  return sum * d / M_PI;
}
}
}
}

#endif  // KERNEL_TESTING_TOOLBOX

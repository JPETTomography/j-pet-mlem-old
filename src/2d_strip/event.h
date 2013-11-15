#pragma once

#include <cmath>

template <typename F> struct event {
  event(F z_u_a, F z_d_a, F dl_a) : z_u(z_u_a), z_d(z_d_a), dl(dl_a) {}
  event() {};
  F z_u;
  F z_d;
  F dl;
};

template <typename F> struct ImageSpaceEventAngle {
  ImageSpaceEventAngle(F y_a, F z_a, F angle_a)
      : y(y_a), z(z_a), angle(angle_a) {};
  const F y;
  const F z;
  const F angle;
};

template <typename F> struct ImageSpaceEventTan {
  ImageSpaceEventTan(F y_a, F z_a, F tan_a) : y(y_a), z(z_a), tan(tan_a) {};
  const F y;
  const F z;
  const F tan;
};

template <typename F>
ImageSpaceEventTan<F> to_tan(const ImageSpaceEventAngle<F>& ea) {
  return ImageSpaceEventTan<F>(ea.y, ea.z, tan(ea.angle));
}

template <typename F>
ImageSpaceEventAngle<F> to_angle(const ImageSpaceEventTan<F>& ea) {
  return ImageSpaceEventAngle<F>(ea.y, ea.z, atan(ea.tan));
}

template <typename T> struct ellipse_parameters {

  T x, y, a, b;
  T angle;
  T iter;
};

template <typename T> T event_tan(T z_u, T z_d, T R) {
  return (z_u - z_d) / (T(2.0) * R);
}
template <typename T> T event_y(T dl, T tan_event) {
  return -T(0.5) * (dl / std::sqrt(T(1.0) + (tan_event * tan_event)));
}
template <typename T> T event_z(T z_u, T z_d, T y, T tan_event) {
  return T(0.5) * (z_u + z_d + (T(2.0) * y * tan_event));
}

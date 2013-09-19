#pragma once

template <typename T> struct event {

  T z_u;
  T z_d;
  T dl;
};

template <typename T> struct ellipse_parameters {

  T x, y, a, b;
  T angle;
  T iter;
};

template <typename T> T event_tan(T z_u, T z_d, T R) {
  return (z_u - z_d) / (T(2.0) * R);
}
template <typename T> T event_y(T dl, T tan_event) {
  return -T(0.5) * (dl / sqrt(T(1.0) + (tan_event * tan_event)));
}
template <typename T> T event_z(T z_u, T z_d, T y, T tan_event) {
  return T(0.5) * (z_u + z_d + (T(2.0) * y * tan_event));
}

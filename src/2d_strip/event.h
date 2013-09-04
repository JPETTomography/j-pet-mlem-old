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

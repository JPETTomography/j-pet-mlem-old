#pragma once

#define DEBUG 0
#define DEBUG_KERNEL 0
#define DEBUG_BBOX 1
#define DEBUG_PIXEL_IN_ELLIPSE 0

template <typename T> struct event {

  T z_u;
  T z_d;
  T dl;
};

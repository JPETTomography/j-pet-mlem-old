#pragma once

#define DEBUG 0
#define DEBUG_KERNEL 1

template <typename T> struct event {

  T z_u;
  T z_d;
  T dl;
};

#pragma once

#include "../point.h"
#include "../detector.h"

template <typename F> class ToR {
 public:
  ToR(Point<F> c1, F angle1, F w1, F h1, Point<F> c2, F angle2, F w2, F h2)
      : d_({ { h1, w1 }, { h2, w2 } }),
        c_{ c1, c2 },
        angle_{ angle1, angle2 },
        w_{ w1, w2 },
        h_{ h1, h2 } {
    for (int i = 0; i < 2; i++) {
      d_[i].rotated(angle_[i]);
      d_[i] += c_[i];
    }
  }

  Point<F> center(int i) const { return c_[i]; }
  F width(int i) const { return w_[i]; }
  F height(int i) const { return h_[i]; }
  F angle(int i) const { return angle_[i]; }
  Detector<F> detector(int i) const { return d_[i]; }

 private:
  Detector<F> d_[2];
  Point<F> c_[2];
  F angle_[2];
  F w_[2], h_[2];
};

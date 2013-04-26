#pragma once

#include"../point.h"
#include"../detector.h"

template <typename F> class ToR {
 public:
  ToR(Point<F> c1, F angle1, F w1, F h1, Point<F> c2, F angle2, F w2, F h2)
      : d_({
    { h1, w1 }
    , { h2, w2 }
  }),
        c_ {
    c1, c2
  }
  , w_ { w1, w2 }
  , h_ { h1, h2 }
  {
    d_[0].rotated(angle1);
    d_[0] += c1;
    d_[1].rotated(angle2);
    d_[1] += c2;
  }
  ;

 private:
  Detector<F> d_[2];
  Point<F> c_[2];
  F w_[2], h_[2];

};

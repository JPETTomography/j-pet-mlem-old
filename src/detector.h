#pragma once

#include <utility>

#include "event.h"

template<typename S, typename F> class Ring2DDetector {
public:
  Ring2DDetector(S n_detectors, S n_pixels):
    n_detectors_(n_detectors),
    n_pixels_x_(n_pixels),
    n_pixels_y_(n_pixels),
    pixel_size_x_(2.0/n_pixels),
    pixel_size_y_(2.0/n_pixels),
    radius_(sqrt(2.0)) {}

  S n_detectors() const {return n_detectors_; }
  S n_pixels_x() const {return n_pixels_x_; }
  S n_pixels_y() const {return n_pixels_y_; }
  F pixel_size_x() const {return pixel_size_x_; }
  F pixel_size_y() const {return pixel_size_y_; }
  F radius() const {return radius_; }

private:
  S n_detectors_;
  S n_pixels_x_;
  S n_pixels_y_;
  F pixel_size_x_;
  F pixel_size_y_;
  F radius_;
};

template<typename F>
std::pair<F, F> tof(const event<F> ev, F R2) {
  F s, c;
  sincos(ev.phi(), s, c);
  F x = ev.x();
  F y = ev.y();

  F b = x*c +y*s;

  F delta = sqrt(b*b+R2-x*x-y*y);

  F t1=-b+delta;
  F t2=-b-delta;

  return std::make_pair(t1, t2);
};

template<typename T, typename F>
struct Lor {
  Lor() {};
  Lor(T f, T s, F p = 0.0):first(f), second(s), count(p) {
    if(s<f) std::swap(s, f);
  }
  T first;
  T second;
  F count;
};

template<typename F>
std::pair<short, short>   lor(const std::pair<F, F> &time,
           const event<F> ev,
           F R2, int n) {

  F a_step = 2.0*M_PI/n;

  F s, c;
  sincos(ev.phi(), s, c);
  F x = ev.x();
  F y = ev.y();

  F x_hit = x+time.first*c;
  F y_hit = y+time.first*s;

  F phi_hit= atan2(y_hit, x_hit);
  short l1=(short)floor(phi_hit/a_step);

  if(l1<0)
    l1+=n;

  x_hit = x+time.second*c;
  y_hit = y+time.second*s;

  phi_hit= atan2(y_hit, x_hit);
  short l2=(short)floor(phi_hit/a_step);

  if(l2<0)
    l2+=n;

  if(l1<=l2)
    return std::make_pair(l1, l2);
  else
    return std::make_pair(l2, l1);
};

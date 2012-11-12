#pragma once

#include "taus.h"

#include "tof_event.h"
#include "tof_detector.h"
#include "phantom.h"

template<typename F, typename D> class ToF_Monte_Carlo  {
public:
  ToF_Monte_Carlo(const D &detector, int n_gen = 1):
    detector_(detector), taus_(n_gen) {};

  void gen_seeds(unsigned long seed) {taus_.gen_seeds(seed);}

  void sc(double x, double *s, double *c) { sincos(x, s, c); }
  void sc(float x, float *s, float *c)    { sincosf(x, s, c); }

  ToF_Track_2D<F> add_noise(const  ToF_Track_2D<F> &track, int gen = 0) {
    F r1 = taus_.rnd(gen);
    F r2 = taus_.rnd(gen);
    F angle = 2.0*M_PI*r1;
    F l = sqrt(-2.0*log(r2));
    F s, c;
    sc(angle, &s, &c);
    F g1 = l*s;
    F g2 = l*c;
    F z_up = track.z_up()+detector_.sigma_x()*g1;
    F z_dn = track.z_dn()+detector_.sigma_x()*g2;
    r1 = taus_.rnd(gen);
    r2 = taus_.rnd(gen);
    angle = 2.0*M_PI*r1;
    l = sqrt(-2.0*log(r2));
    sc(angle, &s, &c);
    g1 = l*s;

    F dl = track.dl()+detector_.sigma_l()*g1;

    return ToF_Track_2D<F>(z_up, z_dn, dl);
  }
  ToF_Event_2D<F> add_noise(const  ToF_Event_2D<F> &event, int gen = 0) {
    return detector_.fromPS(add_noise(detector_.toPS(event)));
  }

  template<typename In, typename Out>
  int  add_noise(In first, In last, Out out, int gen = 0) {
    int count = 0;
    while(first != last) {
      *out = add_noise(*first);
      ++out;
      ++first;
      ++count;
    }
    return count;
  }

  template<typename In, typename Out>
  int  add_noise_to_detected(In first, In last, Out out, int gen = 0) {
    int count = 0;
    while(first != last) {
      ToF_Track_2D<F> track = detector_.toPS(*first);
      if(detector_.detected(track)) {
      *out = add_noise(*first);
      ++out;
      ++count;
      }
      ++first;
    }
    return count;
  }

  ToF_Event_2D<F> emit_from_point(F x, F y, int gen = 0)  {
    F phi = 2.0*M_PI*taus_.rnd(gen);
    return  ToF_Event_2D<F>(x, y, tan(phi));
  }

  int
  fill_with_events_from_single_point(std::vector< ToF_Event_2D<F> > &events, F x, F y, int n) {
    for(int i = 0;i<n;++i) {
      events[i] = emit_from_point(x, y);
    }
    return n;
  }

  template<typename Out> int
  fill_with_events_from_single_point(Out out, F x, F y, int n, int gen = 0) {
    for(int i = 0;i<n;++i, ++out) {
      *out = emit_from_point(x, y, gen);
    }
    return n;

  }

  template<typename Out> int
  fill_with_detected_events_from_single_point(Out out, F x, F y, int n, int gen = 0) {
    int count = 0;
    for(int i = 0;i<n;++i) {
      ToF_Event_2D<F> event = emit_from_point(x, y, gen);
      if(detector_.detected(event)) {
    *out = event;
    out++;
    count++;
  }
    }
    return count;
  }
  template<typename Out> int
  fill_with_events_from_phantom(Out out, Phantom *phantom, F ll_x, F ll_y, F ur_x, F ur_y, int n, int gen = 0 ) {
    F dx = ur_x-ll_x;
    F dy = ur_y-ll_y;
    int count = 0;
    for(int i = 0;i<n;++i) {
      F x = ll_x+dx*taus_.rnd(gen);
      F y = ll_y+dy*taus_.rnd(gen);
      F rnd = taus_.rnd(gen);
      if(phantom->emit(x, y, rnd)) {
  *out = emit_from_point(x, y, gen);
    ++count;
    ++out;
      }
    }
    return count;
  }
  private:
    taus_array taus_;
    D detector_;

  };

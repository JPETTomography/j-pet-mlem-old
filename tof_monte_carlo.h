#ifndef __TOF_MONTE_CARLO_H__
#define __TOF_MONTE_CARLO_H__

#include"taus.h"

#include"tof_event.h"
#include"tof_detector.h"
#include"phantom.h"


template<typename F,  typename  D> class ToF_Monte_Carlo  {
public:
  ToF_Monte_Carlo(const D &detector, int n_gen=1):
    detector_(detector),
    taus_(n_gen) {};
  
  void gen_seeds(unsigned long seed) {taus_.gen_seeds(seed);}

  ToF_Track_2D<F> add_noise(const  ToF_Track_2D<F> &track,int gen = 0) {
    F r1=taus_.rnd(gen);
    F r2=taus_.rnd(gen);
    F angle=2.0*M_PI*r1;
    F l=sqrt(-2.0*log(r2));
    F s,c;
    sincos(angle,&s,&c);
    F g1=l*s;
    F g2=l*c;

    std::cerr<<"gauss  "<<g1<<std::endl;

    F z_up=track.z_up()+detector_.sigma_x()*g1;
    F z_dn=track.z_dn()+detector_.sigma_x()*g2;
    r1=taus_.rnd(gen);
    r2=taus_.rnd(gen);
    angle=2.0*M_PI*r1;
    l=sqrt(-2.0*log(r2));
    s,c;
    sincos(angle,&s,&c);
    g1=l*s;
    
    F dl=track.dl()+detector_.sigma_l()*g1;

    return ToF_Track_2D<F>(z_up,z_dn,dl);
  }
  

  ToF_Event_2D<F> add_noise(const  ToF_Event_2D<F> &event,int gen = 0) {
    return detector_.fromPS(add_noise(detector_.toPS(event)));
  }

private:
  taus_array taus_;
  D detector_;

};

#endif

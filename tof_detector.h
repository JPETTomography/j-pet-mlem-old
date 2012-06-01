#ifndef __TOF_DETECTOR_H__
#define __TOF_DETECTOR_H__

#include"tausworthe.h"
#include"tof_event.h"

template<typename F>
class ToF_Detector_2D {
public:

  ToF_Detector_2D(F R, F L):R_(R),L_(L) {};
  

  F sigma_x() const {return sigma_x_;}
  F sigma_l() const {return sigma_l_;}

  void set_sigma(F sigma_x, F sigma_l) {
    sigma_x_=sigma_x;
    sigma_l_=sigma_l;
  }

  ToF_Event_2D<F> fromPS(const ToF_Track_2D<F> &track) {
    double tan=(track.z_up()-track.z_dn())/((F)2.0*R_);
    double y=-0.5*track.dl()/sqrt(tan*tan+1);
    double z=(F)(0.5)*(track.z_up()+track.z_dn())+y*tan;

    return ToF_Event_2D<F>(z,y,tan);
  }

  ToF_Track_2D<F> toPS(const ToF_Event_2D<F> &event) {
    double z_up=event.z()+(R_-event.y())*event.tan();
    double z_dn=event.z()-(R_+event.y())*event.tan();
    double dl=-2.0*event.y()*sqrt(event.tan()*event.tan()+1);
    return ToF_Track_2D<F>(z_up,z_dn,dl);
}
 
   

 
  typedef F float_t;
  
private:
  F R_;
  F L_;
  F sigma_x_;
  F sigma_l_;

};

#endif

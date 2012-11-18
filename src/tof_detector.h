#ifndef __TOF_DETECTOR_H__
#define __TOF_DETECTOR_H__

#include "random.h"
#include "tof_event.h"

template<typename F>
class ToF_Detector_2D {
public:
  typedef F float_t;
  typedef ToF_Track_2D<F> track_t;
  typedef ToF_Event_2D<F> event_t;
 ToF_Detector_2D(F R, F L):R_(R), L_(L), L_h_(L/2.0) {};

  F R() const {return R_;}
  F L() const {return L_;}

  F sigma_x() const {return sigma_x_;}
  F sigma_l() const {return sigma_l_;}

  void set_sigma(F sigma_x, F sigma_l) {
    sigma_x_ = sigma_x;
    sigma_l_ = sigma_l;
  }

  ToF_Event_2D<F> fromPS(const ToF_Track_2D<F> &track) {
    double tan = (track.z_up()-track.z_dn())/((F)2.0*R_);
    double y = -0.5*track.dl()/sqrt(tan*tan+1);
    double z = (F)(0.5)*(track.z_up()+track.z_dn())+y*tan;

    return ToF_Event_2D<F>(z, y, tan);
  }

  ToF_Track_2D<F> toPS(const ToF_Event_2D<F> &event) {
    double z_up = event.z()+(R_-event.y())*event.tan();
    double z_dn = event.z()-(R_+event.y())*event.tan();
    double dl = -2.0*event.y()*sqrt(event.tan()*event.tan()+1);
    return ToF_Track_2D<F>(z_up, z_dn, dl);
}

  bool detected(const ToF_Track_2D<F> &track) {
    F z_up = track.z_up();
    F z_dn = track.z_dn();

    return z_up < L_h_ && z_up > -L_h_ && z_dn < L_h_ && z_dn>-L_h_ ;

  }
  bool detected(const ToF_Event_2D<F> &event) {
    track_t track = toPS(event);
    F z_up = track.z_up();
    F z_dn = track.z_dn();

    return z_up < L_h_ && z_up > -L_h_ && z_dn < L_h_ && z_dn>-L_h_ ;

  }
  F acceptance(F x, F y) const {
    F z_r = ((R_+y)*L_h_+2.0*x*R_)/(R_-y);
    F phi_l, phi_r;
    if(z_r<L_h_)
      phi_l = atan( (L_h_+x)/(R_-y) );
    else
      phi_l = atan( (L_h_-x)/(R_+y) );

    F z_l = (-(R_+y)*L_h_+2.0*x*R_)/(R_-y);

    if(z_l>-L_h_)
      phi_r = atan( (L_h_-x)/(R_-y) );
    else
      phi_r = atan( (L_h_+x)/(R_+y) );
    return (phi_r+phi_l)/M_PI;
  }
private:
  F R_;
  F L_;
  F L_h_;
  F sigma_x_;
  F sigma_l_;
};

#endif

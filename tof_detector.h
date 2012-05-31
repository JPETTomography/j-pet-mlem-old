#ifndef __TOF_DETECTOR_H__
#define __TOF_DETECTOR_H__

#include"tof_event.h"

template<typename F>
class ToF_Detector_2D {
public:

  ToF_Detector_2D(F R, F L):R_(R),L_(L) {};

  ToF_Event_2D<F> fromPS(F z_up, F z_dn, F dl) {
    double tan=(z_up-z_dn)/((F)2.0*R_);
    double z=(F)(0.5)*(z_up+z_dn);
    double y=-dl/sqrt(tan*tan+1);
    return ToF_Event_2D<F>(z,y,tan);
  }

  ToF_Track_2D<F> toPS(F z_up, F z_dn, F dl,F R) {
    double tan=(z_up-z_dn)/((F)2.0*R);
    double z=(F)(0.5)*(z_up+z_dn);
    double y=-dl/sqrt(event.tan_*event.tan_+1);
    return ToF_Event_2D(z,y,tan);
  }
  
  typedef F float_t;
  
private:
  F R_;
  F L_;


};

#endif

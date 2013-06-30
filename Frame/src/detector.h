#pragma once
#include<cmath>

const double Degree = M_PI/180.0;

template<typename F> 
class EventImageAngle {
public:
  EventImageAngle(F y, F z, F theta): y_(y), z_(z), theta_(theta) {};

  F  y()     const {return y_;}
  F  z()     const {return z_;}
  F  theta() const {return theta_;}

private:
  const F y_, z_, theta_;
};

template<typename F> 
class EventImageTan {
public:
  EventImageTan(F y, F z, F tangent_): y_(y), z_(z), tangent_(tangent_) {};

  F  y()     const {return y_;}
  F  z()     const {return z_;}
  F  tan() const {return tangent_;}

private:
  const F y_, z_, tangent_;
};

template<typename F>
EventImageAngle<F> EventImageTanToAngle(EventImageTan<F> event) {
  return EventImageAngle<F>(event.y(), event.z(), ::atan( event.tan() ) );
}

template<typename F>
EventImageTan<F> EventImageAngleToTan(EventImageAngle<F> event) {
  return EventImageTan<F>(event.y(), event.z(), ::tan( event.theta() ) );
}


template<typename F>
class EventDetector {
public:
  EventDetector(F z_up, F z_dn, F dl): z_up_(z_up), z_dn_(z_dn), dl_(dl) {};

  F  z_up()   const {return z_up_;}
  F  z_dn()   const {return z_dn_;}
  F  dl()     const {return dl_;}

private:
  const F z_up_, z_dn_, dl_;
};



template<typename F> 
class Detector {
public:
  Detector(F R, F L):R_(R), L_(L) {};
  
  F R() const  {return R_;}
  F L() const  {return L_;}

  EventDetector<F> EventImageTanToDetector(const EventImageTan<F> &event) {
    F z_up = event.z() + (R()-event.y())*event.tan();
    F z_dn = event.z() - (R()+event.y())*event.tan();
    F dl   = -2.0*event.y()*::sqrt(1+event.tan()*event.tan());
    return EventDetector<F>(z_up, z_dn, dl);
  }; 

 EventDetector<F> EventImageAngleToDetector(const EventImageAngle<F> &event) {
   return EventImageTanToDetector(EventImageAngleToTan(event));
  }; 



  EventImageTan<F> EventDetectorToImageTan(const EventDetector<F> &event) {
    F tan = ( event.z_up()-event.z_dn() )/(2.0 * R() );
    F y   = -0.5*event.dl()/::sqrt(1+tan*tan);
    F z   = 0.5*(event.z_up()+event.z_dn())+y*tan;
    return EventImageTan<F>(y, z, tan);
  }

  EventImageAngle<F> EventDetectorToImageAngle(const EventDetector<F> &event) {
    return EventImageTanToAngle(EventDetectorToImageTan(event));
  }

private:
  const F R_,L_;
} ;

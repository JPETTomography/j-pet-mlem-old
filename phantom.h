#ifndef __PHANTOM_H_
#define __PHANTOM_H__

#include<cmath>

class ElipticalRegion {
 public:

 ElipticalRegion(double x, double y, double a, double b, double phi,double act):
  x_(x),y_(y),a_(a),b_(b),phi_(phi),activity_(act) {
    sincos(&sin_,&cos_);
    inv_a2_=1.0/(a_*a_);
    inv_b2_=1.0/(b_*b_);
  };
  double activity() const {return activity_;}
  bool in(double x,double y) const {

    double dx=x-x_;
    double dy=y-y_;

    double r_x= dx * cos_+dy*sin_;
    double r_y=-dx*sin_ +dy*cos_;

    double r2=r_x*r_x*inv_a2_+r_y*r_y*inv_b2_;
    
    return r2<1.0;

  }

  bool emited(double r) const {
    return activity() > r;
  }


 private:
  double x_;
  double y_;
  double a_;
  double b_;
  
  double phi_;


  double activity_;

  double inv_a2_;
  double inv_b2_;
  
  double sin_;
  double cos_;

}

#endif

#ifndef __PHANTOM_H__
#define __PHANTOM_H__

#include<cmath>

#include<vector>

class EllipticalRegion {
 public:

 EllipticalRegion(double x, double y, double a, double b, double phi,double act):
  x_(x),y_(y),a_(a),b_(b),phi_(phi),activity_(act) {
    sincos(phi,&sin_,&cos_);
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

};

class Phantom {



 typedef   std::vector<EllipticalRegion *> container;

 public:
 
  void addRegion(double x, double y, double a, double b, double phi,double act) {
    regions_.push_back(new EllipticalRegion(x,y,a,b,phi,act));
  }


  double activity(double x, double y) const {
    container::const_reverse_iterator rit=regions_.rbegin();
    for(;rit!=regions_.rend();++rit) {

      //std::cerr<<(*rit)->activity()<<std::endl;
      if((*rit)->in(x,y)) {

	//std::cerr<<"RETURNING "<<(*rit)->activity()<<std::endl;
	return (*rit)->activity();
      }
    }

    return 0.0;

  }


  bool emit(double x, double y, double rnd) const {
    return activity(x,y)> rnd;
  }


  ~Phantom() {
    container::iterator rit=regions_.begin();
    for(;rit!=regions_.end();++rit) {
      delete (*rit);
    }
  }
  
private:
  container regions_; 

  

};

#endif

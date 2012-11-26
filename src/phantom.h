#pragma once

#include <cmath>
#include <vector>
#include "point.h"

#ifdef __APPLE__
static inline void sincos(double a, double *s, double *c) { *s = sin(a);  *c = cos(a);  };
static inline void sincosf(float a, float *s,  float *c)  { *s = sinf(a); *c = cosf(a); };
#endif

class EllipticalRegion {
 public:

 EllipticalRegion(double x, double y, double a, double b, double phi, double act):
  x_(x), y_(y), a_(a), b_(b), phi_(phi), activity_(act) {
    sincos(phi, &sin_, &cos_);
    inv_a2_ = 1.0/(a_*a_);
    inv_b2_ = 1.0/(b_*b_);
  };

  double activity() const {return activity_;}
  bool in(double x, double y) const {

    double dx = x-x_;
    double dy = y-y_;

    double r_x = dx * cos_+dy*sin_;
    double r_y = -dx*sin_ +dy*cos_;

    double r2 = r_x*r_x*inv_a2_+r_y*r_y*inv_b2_;

    return r2<1.0;

  }
  double x() const {return x_;}
  double y() const {return y_;}
  double a() const {return a_;}
  double b() const {return b_;}
  double phi() const {return phi_;}

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

 typedef container::const_iterator const_iterator;
  Phantom(){};
#if 1
 Phantom(double ll_x, double ll_y, double ur_x, double ur_y):
 ll_x_(ll_x), ll_y_(ll_y), ur_x_(ur_x), ur_y_(ur_y) {}
#endif

  size_t n_regions() const {return regions_.size();}

  void addRegion(double x, double y, double a, double b, double phi, double act) {


    regions_.push_back(new EllipticalRegion(x, y, a, b, phi, act));
  }
  double activity(double x, double y) const {
    container::const_reverse_iterator rit = regions_.rbegin();
    for(;rit != regions_.rend();++rit) {

      if((*rit)->in(x, y)) {
  return (*rit)->activity();
      }
    }

    return 0.0;

  }
  bool emit(double x, double y, double rnd) const {
    return activity(x, y)> rnd;
  }
  const_iterator begin() const {return regions_.begin();}
  const_iterator end() const {return regions_.end();}
  ~Phantom() {
    container::iterator rit = regions_.begin();
    for(;rit != regions_.end();++rit) {
      delete (*rit);
    }
  }

private:
  double ll_x_, ll_y_;
  double ur_x_, ur_y_;

  container regions_;
};


template<typename F>
struct point_source_t {
  point_source_t(F x, F y, F intensity_a): p(x,y),intensity(intensity_a) {}
  point<F> p;
  F intensity;
};


template<typename F>
class  point_sources_t {
public:


  size_t n_sources() const {return sources_.size();}

  void add(F x , F y, F intensity) {
    sources_.push_back(point_source_t<F>(x,y,intensity)); 
  }


  void normalize() {
    total_=0.0;
    for(auto it=sources_.begin();
        it!=sources_.end();
        ++it) {
      total_+=(*it).intensity;
    }

    F cumulant=0.0;
    cumulant_.clear();
    for(auto it=sources_.begin();
        it!=sources_.end();
        ++it) {
      cumulant+=(*it).intensity/=total_;
      cumulant_.push_back(cumulant);
    }
  }


  point<F> draw(F rng) {
    
    int i=0;

    while(rng>cumulant_[i]) {
      ++i;
    }
    return sources_[i].p;
  }

private:
  std::vector<point_source_t<F> > sources_;
  std::vector<F> cumulant_;
  F total_;
};

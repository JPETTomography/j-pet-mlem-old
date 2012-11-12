#pragma once

#include "polygon.h"

template <typename F = double>
class detector : public polygon<F> {
public:
  typedef F angle_type;
  typedef typename polygon<F>::point_type point_type;
  detector(F w, F h)
  {
    this->push_back( {  w/2., h/2.} );
    this->push_back( {  w/2.,-h/2.} );
    this->push_back( { -w/2.,-h/2.} );
    this->push_back( { -w/2., h/2.} );
  }

  detector rotated(angle_type phi) {
    detector r;
    for(auto it = this->begin(); it != this->end(); ++it) {
      r.push_back( (*it).rotated(phi) );
    }
    return r;
  }

  detector & operator += (point_type t) {
    for(auto it = this->begin(); it != this->end(); ++it) *it += t;
    return *this;
  }

private:
  detector() {}
};

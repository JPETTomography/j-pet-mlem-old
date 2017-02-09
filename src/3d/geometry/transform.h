#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"

template <typename F>
PET3D::Vector<F> rotate(const PET3D::Vector<F>& v,
                        F theta,
                        const PET3D::Vector<F>& d) {
  return v;
}

template <typename F>
PET3D::Point<F> rotate(const PET3D::Point<F>& p,
                       F theta,
                       const PET3D::Vector<F>& d,
                       const PET3D::Point<F>& c) {
  return p;
}

#endif  // TRANSFORM_H

#ifndef TRANSFORM_H
#define TRANSFORM_H

#include "3d/geometry/vector.h"
#include "3d/geometry/point.h"

template <typename F>
PET3D::Vector<F> rotate(const PET3D::Vector<F>& v,
                        F theta,
                        const PET3D::Vector<F>& d) {
  F c = std::cos(theta);
  F s = std::sin(theta);
  F R[3][3] = { { c + d.x * d.x * (1 - c),
                  d.x * d.y * (1 - c) - d.z * s,
                  d.x * d.z * (1 - c) + d.y * s },
                { d.x * d.y * (1 - c) + d.z * s,
                  c + d.y * d.y * (1 - c),
                  d.y * d.z * (1 - c) - d.x * s },
                { d.x * d.z * (1 - c) - d.y * s,
                  d.y * d.z * (1 - c) + d.x * s,
                  c + d.z * d.z * (1 - c) } };
  PET3D::Vector<F> w(0, 0, 0);

  w.x = R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z;
  w.y = R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z;
  w.z = R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z;

  return w;
}

template <typename F>
PET3D::Point<F> rotate(const PET3D::Point<F>& p,
                       F theta,
                       const PET3D::Vector<F>& d,
                       const PET3D::Point<F>& c) {
  return p;
}

#endif  // TRANSFORM_H

#pragma once

#include <cmath>
#if _MSC_VER
#include <array>
#endif

#include "geometry/point.h"
#include "util/cuda/compat.h"

template <typename FType = double> class Kernel {
 public:
  typedef FType F;
  typedef ::Point<F> Point;
  typedef FType FVec[3];

  const F inv_pow_sigma_z;
  const F inv_pow_sigma_dl;
  const F sqrt_det_cor_mat;

  _ Kernel(F sigma_z, F sigma_dl)
      : inv_pow_sigma_z(1 / (sigma_z * sigma_z)),
        inv_pow_sigma_dl(1 / (sigma_dl * sigma_dl)),
        sqrt_det_cor_mat(inv_pow_sigma_z *  //
                         inv_pow_sigma_z *  //
                         inv_pow_sigma_dl) {}

 private:
  _ F multiply_by_inv_cor_mat(const FVec vec_a, const FVec vec_b) {
    return vec_a[0] * inv_pow_sigma_z * vec_b[0] +
           vec_a[1] * inv_pow_sigma_z * vec_b[1] +
           vec_a[2] * inv_pow_sigma_dl * vec_b[2];
  }

 public:
  _ F test(const F y,
           const F z,
           const Point pixel_center,
           const F dl,
           const F sigma) {

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    return (INV_POW_TWO_PI * (1 / (sigma * dl))) *
           compat::exp(F(-0.5) *
                       (compat::pow((pixel_center.y - y) / dl, F(2)) +
                        compat::pow((pixel_center.x - z) / sigma, F(2))));
  }

  _ F operator()(const F y,
                 const F tan,
                 const F inv_cos,
                 const F pow_inv_cos,  // what power?
                 const F R,
                 const Point pixel_center) {

    FVec vec_o;
    FVec vec_a;
    FVec vec_b;

    vec_o[0] = -(pixel_center.y + y - R) * tan * pow_inv_cos;
    vec_o[1] = -(pixel_center.y + y + R) * tan * pow_inv_cos;
    vec_o[2] = -(pixel_center.y + y) * inv_cos * (1 + 2 * (tan * tan));

    vec_a[0] = -(pixel_center.y + y - R) * pow_inv_cos;
    vec_a[1] = -(pixel_center.y + y + R) * pow_inv_cos;
    vec_a[2] = -2 * (pixel_center.y + y) * (inv_cos * tan);

    vec_b[0] = pixel_center.x - (pixel_center.y * tan);
    vec_b[1] = pixel_center.x - (pixel_center.y * tan);
    vec_b[2] = -2 * pixel_center.y * inv_cos;

    F a_ic_a = multiply_by_inv_cor_mat(vec_a, vec_a);
    F b_ic_a = multiply_by_inv_cor_mat(vec_b, vec_a);
    F b_ic_b = multiply_by_inv_cor_mat(vec_b, vec_b);
    F o_ic_b = multiply_by_inv_cor_mat(vec_o, vec_b);

    F norm = a_ic_a + (2 * o_ic_b);

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    F element_before_exp =
        INV_POW_TWO_PI * (sqrt_det_cor_mat / compat::sqrt(norm));

    F exp_element = -F(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    F exp = compat::exp(exp_element);

    return element_before_exp * exp;
  }

  // TODO: Ellipse bounding box is actually a property of the kernel, not
  //      detector.

  _ void ellipse_bb(F angle,
                    F tan,
                    F& sec,     // out
                    F& sec_sq,  // out
                    F& A,       // out
                    F& B,       // out
                    F& C,       // out
                    F& bb_y,    // out
                    F& bb_z     // out
                    ) const {

    F cos = compat::cos(angle);
    sec = 1 / cos;
    sec_sq = sec * sec;

    A = (4 / (cos * cos)) * inv_pow_sigma_dl + 2 * tan * tan * inv_pow_sigma_z;
    B = -4 * tan * inv_pow_sigma_z;
    C = 2 * inv_pow_sigma_z;
    F B_2 = (B / 2) * (B / 2);

    bb_y = this->bb_y(A, C, B_2);
    bb_z = this->bb_z(A, C, B_2);
  }

  _ bool in_ellipse(F A, F B, F C, Point ellipse_center, Point p) const {

    F dy = p.y - ellipse_center.y;
    F dz = p.x - ellipse_center.x;

    return (A * dy * dy) + (B * dy * dz) + (C * dz * dz) <= 9;
  }

  _ F bb_z(F A, F C, F B_2) const { return 3 / compat::sqrt(C - (B_2 / A)); }
  _ F bb_y(F A, F C, F B_2) const { return 3 / compat::sqrt(A - (B_2 / C)); }
};

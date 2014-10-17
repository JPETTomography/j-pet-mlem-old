#pragma once

#include <cmath>

#include "geometry/point.h"
#include "util/cuda/compat.h"

template <typename FType = double> class Kernel {
 public:
  typedef FType F;
  typedef ::Point<F> Point;
  typedef struct { F p, q, r; } FVec;

 private:
  const F inv_pow_sigma_z;
  const F inv_pow_sigma_dl;
  const F sqrt_det_cor_mat;

 public:
  _ Kernel(F sigma_z, F sigma_dl)
      : inv_pow_sigma_z(1 / (sigma_z * sigma_z)),
        inv_pow_sigma_dl(1 / (sigma_dl * sigma_dl)),
        sqrt_det_cor_mat(inv_pow_sigma_z *  //
                         inv_pow_sigma_z *  //
                         inv_pow_sigma_dl) {}

  _ F test(const F y,
           const F z,
           const Point pixel_center,
           const F dl,
           const F sigma) const {

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    return (INV_POW_TWO_PI * (1 / (sigma * dl))) *
           compat::exp(F(-0.5) *
                       (compat::pow((pixel_center.y - y) / dl, 2) +
                        compat::pow((pixel_center.x - z) / sigma, 2)));
  }

  _ F operator()(const F y,
                 const F cos,
                 const F tan,
                 const F R,
                 const Point pixel_center) {

    F sec = 1 / cos;
    F sec_sq = sec * sec;

    FVec vec_o = { -(pixel_center.y + y - R) * tan * sec_sq,
                   -(pixel_center.y + y + R) * tan * sec_sq,
                   -(pixel_center.y + y) * sec * (1 + 2 * (tan * tan)) };
    FVec vec_a = { -(pixel_center.y + y - R) * sec_sq,
                   -(pixel_center.y + y + R) * sec_sq,
                   -2 * (pixel_center.y + y) * (sec * tan) };
    FVec vec_b = { pixel_center.x - (pixel_center.y * tan),
                   pixel_center.x - (pixel_center.y * tan),
                   -2 * pixel_center.y * sec };

    F a_ic_a = multiply_inv_cor_mat(vec_a, vec_a);
    F b_ic_a = multiply_inv_cor_mat(vec_b, vec_a);
    F b_ic_b = multiply_inv_cor_mat(vec_b, vec_b);
    F o_ic_b = multiply_inv_cor_mat(vec_o, vec_b);

    F norm = a_ic_a + (2 * o_ic_b);

    const F INV_POW_TWO_PI = F(1 / (2 * M_PI * M_PI));

    F element_before_exp =
        INV_POW_TWO_PI * (sqrt_det_cor_mat / compat::sqrt(norm));

    F exp_element = -F(0.5) * (b_ic_b - ((b_ic_a * b_ic_a) / norm));

    F exp = compat::exp(exp_element);

    return element_before_exp * exp;
  }

  _ void ellipse_bb(F cos,
                    F tan,
                    F& A,     // out
                    F& B,     // out
                    F& C,     // out
                    F& bb_y,  // out
                    F& bb_z   // out
                    ) const {

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

 private:
  _ F multiply_inv_cor_mat(const FVec vec_a, const FVec vec_b) const {
    return vec_a.p * inv_pow_sigma_z * vec_b.p +
           vec_a.q * inv_pow_sigma_z * vec_b.q +
           vec_a.r * inv_pow_sigma_dl * vec_b.r;
  }

  _ F bb_z(F A, F C, F B_2) const { return 3 / compat::sqrt(C - (B_2 / A)); }
  _ F bb_y(F A, F C, F B_2) const { return 3 / compat::sqrt(A - (B_2 / C)); }
};

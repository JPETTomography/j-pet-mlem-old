#pragma once

#include <cmath>

#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Strip {

/// Analytic emission probability approximation at given distance to observed
template <typename FType> class GaussianKernel {
 public:
  using F = FType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;

 protected:
  const F inv_pow_sigma_z;
  const F inv_pow_sigma_dl;
  const F inv_pow_two_pi_sqrt_det_cor_mat;
  /// \cond PRIVATE
  typedef struct { F p, q, r; } FVec;
  /// \endcond

 public:
  _ GaussianKernel(F sigma_z, F sigma_dl)
      : inv_pow_sigma_z(1 / (sigma_z * sigma_z)),
        inv_pow_sigma_dl(1 / (sigma_dl * sigma_dl)),
        inv_pow_two_pi_sqrt_det_cor_mat(1 / (2 * M_PI * M_PI) *      //
                                        std::sqrt(inv_pow_sigma_z *  //
                                                  inv_pow_sigma_z *  //
                                                  inv_pow_sigma_dl)) {}

  /// Simplified probability (for testing purposes)
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

  /// Probability of emission at given distance to observed event point
  _ F operator()(const F y,             ///< y position of observed event
                 const F tan,           ///< tangent of event angle
                 const F sec,           ///< \f$ \sqrt(1 + tan * tan) \f$
                 const F R,             ///< radius of frame
                 const Vector distance  ///< distance to emission point
                 ) const {

    F sec_sq = sec * sec;

    F y_position = distance.y + y;
    FVec vec_a = { -(y_position - R) * sec_sq,
                   -(y_position + R) * sec_sq,
                   -2 * y_position * sec * tan };

    F vec_b_pq = distance.x - distance.y * tan;
    FVec vec_b = { vec_b_pq, vec_b_pq, -2 * distance.y * sec };

    F a_ic_a = multiply_inv_cor_mat(vec_a, vec_a);
    F b_ic_b = multiply_inv_cor_mat(vec_b, vec_b);

    F norm = a_ic_a;

    F element_before_exp = inv_pow_two_pi_sqrt_det_cor_mat / compat::sqrt(norm);
    F exp_arg = -F(0.5) * (b_ic_b);
    F exp = compat::exp(exp_arg);

    return element_before_exp * exp;
  }

  /// Probability at given point relative to observed event point
  _ F operator()(const Point event,  ///< observed event point
                 const F tan,        ///< tangent of event angle
                 const F sec,        ///< \f$ \sqrt(1 + tan * tan) \f$
                 const F R,          ///< radius of frame
                 const Point p       ///< point of calculation
                 ) const {
    auto distance = p - event;
    return this->operator()(event.y, tan, sec, R, distance);
  }

  /// Normalized probability at given point relative to observed event point
  _ F normalized(const Point event,  ///< observed event point
                 const F tan,        ///< tangent of event angle
                 const F sec,        ///< \f$ \sqrt(1 + tan * tan) \f$
                 const F R,          ///< radius of frame
                 const F L,          ///< length of frame
                 const Point p       ///< point of calculation
                 ) const {
    auto distance = p - event;
    auto unnormalized = this->operator()(event.y, tan, sec, R, distance);
    auto sens = sensitivity(p, R, L);

    return unnormalized / sens;
  }

  /// Returns \f$ 3 \sigma \f$ ellipse and its bounding box.
  _ void ellipse_bb(F tan,    ///< tangent of event angle
                    F& sec,   ///< [out] \f$ \sqrt(1 + tan * tan) \f$
                    F& A,     ///< [out] \f$ A \$ of ellipse equation
                    F& B,     ///< [out] \f$ B \$ of ellipse equation
                    F& C,     ///< [out] \f$ C \$ of ellipse equation
                    F& bb_y,  ///< [out] ellipse bounding box height
                    F& bb_z   ///< [out] ellipse bounding box weight
                    ) const {

    F tan_sq = tan * tan;
    // NOTE: sqrt(1 + tan*tan) =:= 1 / cos(arctan(tan))
    // sec = 1 / compat::cos(compat::atan(tan));
    sec = compat::sqrt(1 + tan * tan);
    A = 4 * sec * sec * inv_pow_sigma_dl + 2 * tan_sq * inv_pow_sigma_z;
    B = -4 * tan * inv_pow_sigma_z;
    C = 2 * inv_pow_sigma_z;
    F B_2 = (B / 2) * (B / 2);

    bb_y = this->bb_y(A, C, B_2);
    bb_z = this->bb_z(A, C, B_2);
  }

  /// Test if given point lies inside of given ellipse.
  _ bool in_ellipse(F A, F B, F C, Point ellipse_center, Point p) const {

    F dy = p.y - ellipse_center.y;
    F dz = p.x - ellipse_center.x;

    return (A * dy * dy) + (B * dy * dz) + (C * dz * dz) <= 9;
  }

  /// Return frame sensitivity at given point.
  _ static F sensitivity(Point p, F R, F L) {
    auto L2 = L / 2;
    auto theta_max = compat::atan(compat::min((L2 - p.x) / (R - p.y),  //
                                              (L2 + p.x) / (R + p.y)));
    auto theta_min = compat::atan(compat::max(-(L2 + p.x) / (R - p.y),  //
                                              (-L2 + p.x) / (R + p.y)));
    auto sens = theta_max - theta_min;

    return sens / (F)M_PI;
  }

 protected:
  _ F multiply_inv_cor_mat(const FVec vec_a, const FVec vec_b) const {
    return vec_a.p * inv_pow_sigma_z * vec_b.p +
           vec_a.q * inv_pow_sigma_z * vec_b.q +
           vec_a.r * inv_pow_sigma_dl * vec_b.r;
  }

 public:  // needed public for tests
  _ F bb_z(F A, F C, F B_2) const { return 3 / compat::sqrt(C - (B_2 / A)); }
  _ F bb_y(F A, F C, F B_2) const { return 3 / compat::sqrt(A - (B_2 / C)); }
};
}  // Strip
}  // PET2D

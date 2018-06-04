#pragma once

#include <cmath>

#include "2d/geometry/point.h"
#include "util/cuda/compat.h"
#include "gaussian_kernel.h"

namespace PET2D {
namespace Strip {

/// Better analytic approximation than GaussianKernel
////
/// This is better aproximation than GaussianKernel, but note that ellipse and
/// bounding box calculation still uses gaussian approximation.
template <typename FType> class AnalyticKernel : public GaussianKernel<FType> {
 public:
  using F = FType;
  using Gaussian = GaussianKernel<F>;
  using Point = typename Gaussian::Point;
  using Vector = typename Gaussian::Vector;
  using FVec = typename Gaussian::FVec;

 public:
  _ AnalyticKernel(F sigma_z, F sigma_dl) : Gaussian(sigma_z, sigma_dl) {}

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
    FVec vec_o = { vec_a.p * tan,  //
                   vec_a.q * tan,  //
                   -y_position * sec * (1 + 2 * tan * tan) };
    F vec_b_pq = distance.x - distance.y * tan;
    FVec vec_b = { vec_b_pq,  //
                   vec_b_pq,  //
                   -2 * distance.y * sec };

    F a_ic_a = this->multiply_inv_cor_mat(vec_a, vec_a);
    F b_ic_a = this->multiply_inv_cor_mat(vec_b, vec_a);
    F b_ic_b = this->multiply_inv_cor_mat(vec_b, vec_b);
    F o_ic_b = this->multiply_inv_cor_mat(vec_o, vec_b);

    F norm = a_ic_a + 2 * o_ic_b;

    F element_before_exp =
        this->inv_pow_two_pi_sqrt_det_cor_mat / compat::sqrt(norm);
    F exp_arg = -F(0.5) * (b_ic_b - b_ic_a * b_ic_a / norm);
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
    return this->operator()(event.y, tan, sec, R, distance) /
           sensitivity(p, R, L);
  }
};
}  // Strip
}  // PET2D

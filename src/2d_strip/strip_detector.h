#pragma once

#if _MSC_VER
#include <array>
#endif

#include "event.h"

#include "geometry/pixel.h"
#include "geometry/point.h"
#include "util/cuda/compat.h"

/// Class responsible for the Strip detector together with the pixel grid
/// inside.
template <typename FType = double> class StripDetector {
 public:
  typedef FType F;
#if !_MSC_VER || __CUDACC__
  typedef FType FVec[3];
#else
  typedef std::array<F, 3> FVec;
#endif
  typedef ::Pixel<> Pixel;
  typedef ::Point<F> Point;

  StripDetector(F radius,
                F scintilator_length,
                int n_y_pixels,
                int n_z_pixels,
                F pixel_height,  // y direction
                F pixel_width,   // z direction
                F sigma_z,
                F sigma_dl,
                F grid_center_y = 0,
                F grid_center_z = 0)
      : radius(radius),
        scintilator_length(scintilator_length),
        n_y_pixels(n_y_pixels),
        n_z_pixels(n_z_pixels),
        total_n_pixels(n_y_pixels * n_z_pixels),
        pixel_width(pixel_width),
        pixel_height(pixel_height),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl),
        grid_center_y(grid_center_y),
        grid_center_z(grid_center_z),
        grid_size_y(n_y_pixels * pixel_height),
        grid_size_z(n_z_pixels * pixel_width),
        grid_ul_y(grid_center_y + F(0.5) * grid_size_y),
        grid_ul_z(grid_center_z - F(0.5) * grid_size_z),
        inv_pow_sigma_z(1 / (sigma_z * sigma_z)),
        inv_pow_sigma_dl(1 / (sigma_dl * sigma_dl)),
        grid_ul_y_half_h(grid_center_y + F(0.5) * grid_size_y -
                         F(0.5) * pixel_height),
        grid_ul_z_half_w(grid_center_z - F(0.5) * grid_size_z +
                         F(0.5) * pixel_width),
#if !__CUDACC__ && !_MSC_VER
        inv_cor_mat_diag{ 1 / (sigma_z * sigma_z),
                          1 / (sigma_z * sigma_z),
                          1 / (sigma_dl * sigma_dl) },
#elif _MSC_VER
        inv_cor_mat_diag{ { 1 / (sigma_z * sigma_z),
                            1 / (sigma_z * sigma_z),
                            1 / (sigma_dl * sigma_dl) } },
#endif
        half_scintilator_length_(F(0.5) * scintilator_length),
        half_pixel_width_(F(0.5) * pixel_width),
        half_pixel_height_(F(0.5) * pixel_height) {
#if __CUDACC__
    inv_cor_mat_diag[0] = inv_pow_sigma_z;
    inv_cor_mat_diag[1] = inv_pow_sigma_z;
    inv_cor_mat_diag[2] = inv_pow_sigma_dl;
#endif
  }

  Event<F> to_projection_space_tan(const ImageSpaceEventTan<F>& is_event) {
    F z_u = is_event.z + (radius - is_event.y) * is_event.tan;
    F z_d = is_event.z - (radius + is_event.y) * is_event.tan;
    F dl = -2 * is_event.y * sqrt(is_event.tan * is_event.tan + 1);
    return Event<F>(z_u, z_d, dl);
  }

  Event<F> to_projection_space_angle(const ImageSpaceEventAngle<F>& is_ea) {
    return to_angle(to_projection_space_angle(is_ea));
  }

  ImageSpaceEventTan<F> from_projection_space_tan(const Event<F>& event) {
    F tan, y, z;
    event.transform(radius, tan, y, z);
    return ImageSpaceEventTan<F>(y, z, tan);
  }

  ImageSpaceEventAngle<F> from_projection_space_angle(const Event<F>& event) {
    return to_angle(from_projection_space_tan(event));
  }

  _ Point pixel_center(Pixel p) const {
    return Point(grid_ul_z_half_w + p.x * pixel_width,
                 grid_ul_y_half_h - p.y * pixel_height);
  }

  _ Pixel pixel_location(Point p) const {
    return Pixel((grid_ul_y - p.y) / pixel_height,
                 (p.x - grid_ul_z) / pixel_width);
  }

  _ F sensitivity(Point p) const {
    F L_plus = half_scintilator_length_ + p.x;
    F L_minus = half_scintilator_length_ - p.x;
    F R_plus = radius + p.y;
    F R_minus = radius - p.y;

    return F(M_1_PI) *
           (compat::atan(compat::min(L_minus / R_minus, L_plus / R_plus)) -
            compat::atan(compat::max(-L_plus / R_minus, -L_minus / R_plus)));
  }

  _ F sqrt_det_cor_mat() const {
    return compat::sqrt(inv_cor_mat_diag[0] *  //
                        inv_cor_mat_diag[1] *  //
                        inv_cor_mat_diag[2]);
  }

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

    return (A * (dy * dy)) + (B * dy * dz) + (C * (dz * dz)) <= 9;
  }

  const F radius;
  const F scintilator_length;
  const int n_y_pixels;
  const int n_z_pixels;
  const int total_n_pixels;
  const F pixel_width;
  const F pixel_height;
  const F sigma_z;
  const F sigma_dl;
  const F grid_center_y;
  const F grid_center_z;
  const F grid_size_y;
  const F grid_size_z;
  const F grid_ul_y;
  const F grid_ul_z;
  const F inv_pow_sigma_z;
  const F inv_pow_sigma_dl;
  const F grid_ul_y_half_h;
  const F grid_ul_z_half_w;
  const FVec inv_cor_mat_diag;

 private:
  const F half_scintilator_length_;
  const F half_pixel_width_;
  const F half_pixel_height_;

  _ F bb_z(F A, F C, F B_2) const { return 3 / compat::sqrt(C - (B_2 / A)); }
  _ F bb_y(F A, F C, F B_2) const { return 3 / compat::sqrt(A - (B_2 / C)); }
};

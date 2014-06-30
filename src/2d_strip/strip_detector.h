#pragma once
#include "event.h"
/**
 * Class  responsible for the Strip detector together with
 * the pixel grid inside.
 */

template <typename F> class StripDetector {
  static constexpr const F INVERSE_PI = F(M_1_PI);

 public:
  typedef std::pair<int, int> Pixel;
  typedef std::pair<F, F> Point;

  StripDetector(F radius,
                F scintilator_length,
                int n_y_pixels,
                int n_z_pixels,
                F pixel_height,  // y direction
                F pixel_width,   // z direction
                F sigma_z,
                F sigma_dl,
                F grid_center_y = 0.0,
                F grid_center_z = 0.0)
      : radius(radius),
        scintilator_length(scintilator_length),
        n_y_pixels(n_y_pixels),
        n_z_pixels(n_z_pixels),
        pixel_width(pixel_width),
        pixel_height(pixel_height),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl),
        grid_center_y(grid_center_y),
        grid_center_z(grid_center_z),
        grid_size_y_(n_y_pixels * pixel_height),
        grid_size_z_(n_z_pixels * pixel_width),
        grid_ul_y_(grid_center_y + 0.5 * grid_size_y_),
        grid_ul_z_(grid_center_z - 0.5 * grid_size_z_) {
    F s_z_sq = sigma_z * sigma_z;
    F s_dl_sq = sigma_dl * sigma_dl;

    F inv_s_z_sq = 1.0 / s_z_sq;
    F inv_s_dl_sq = 1.0 / s_dl_sq;
    inverse_correlation_matrix_[0][0] = inv_s_z_sq;
    inverse_correlation_matrix_[0][1] = F(0.0f);
    inverse_correlation_matrix_[0][2] = F(0.0f);

    inverse_correlation_matrix_[1][0] = F(0.0f);
    inverse_correlation_matrix_[1][1] = inv_s_z_sq;
    inverse_correlation_matrix_[1][2] = F(0.0f);

    inverse_correlation_matrix_[2][0] = F(0.0f);
    inverse_correlation_matrix_[2][1] = F(0.0f);
    inverse_correlation_matrix_[2][2] = inv_s_dl_sq;
  }

  F half_scintilator_length() const { return F(0.5) * scintilator_length; }

  F inv_c(int i, int j) const { return inverse_correlation_matrix_[i][j]; }

  Event<F> to_projection_space_tan(const ImageSpaceEventTan<F>& ev) {
    F z_u = ev.z + (radius - ev.y) * ev.tan;
    F z_d = ev.z - (radius + ev.y) * ev.tan;
    F dl = -F(2.0) * ev.y * sqrt(ev.tan * ev.tan + F(1.0));
    return Event<F>(z_u, z_d, dl);
  }

  Event<F> to_projection_space_angle(const ImageSpaceEventAngle<F>& img_event) {
    return to_angle(to_projection_space_angle(img_event));
  }

  ImageSpaceEventTan<F> from_projection_space_tan(const Event<F>& ev) {
    F t = event_tan(ev.z_u, ev.z_d, radius);
    F y = event_y(ev.dl, t);
    F z = event_z(ev.z_u, ev.z_d, y, t);
    return ImageSpaceEventTan<F>(y, z, t);
  }

  ImageSpaceEventAngle<F> from_projection_space_angle(const Event<F>& ev) {
    return to_angle(from_projection_space_tan(ev));
  }

  Point pixel_center(int i, int j) {
    return std::make_pair<F>(grid_ul_y_ - i * pixel_height - 0.5 * pixel_height,
                             grid_ul_z_ + j * pixel_width + 0.5 * pixel_width);
  }

  Point pixel_center(Pixel pix) { return pixel_center(pix.first, pix.second); }

  Pixel pixel_location(F y, F z) {
    return std::make_pair<int>(floor((grid_ul_y_ - y) / pixel_height),
                               floor((z - grid_ul_z_) / pixel_width));
  }

  Pixel pixel_location(Point p) { return pixel_location(p.first, p.second); }

  F sensitivity(F y, F z) {
    F L_plus = (half_scintilator_length() + z);
    F L_minus = (half_scintilator_length() - z);
    F R_plus = radius + y;
    F R_minus = radius - y;

    return INVERSE_PI *
           (std::atan(std::min(L_minus / R_minus, L_plus / R_plus)) -
            std::atan(std::max(-L_plus / R_minus, -L_minus / R_plus)));
  }

  const F radius;
  const F scintilator_length;
  const int n_y_pixels;
  const int n_z_pixels;
  const F pixel_width;
  const F pixel_height;
  const F sigma_z;
  const F sigma_dl;
  const F grid_center_y;
  const F grid_center_z;

 public:
  F grid_size_y_;
  F grid_size_z_;
  F grid_ul_y_;
  F grid_ul_z_;

  F inverse_correlation_matrix_[3][3];
};

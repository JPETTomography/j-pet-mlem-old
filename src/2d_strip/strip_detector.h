#pragma once
#include "event.h"
/**
 * Class  responsible for the Strip detector together with
 * the pixel grid inside.
 */

template <typename F> class StripDetector {
  static constexpr const F INVERSE_PI = F(1.0 / M_PI);

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
      : radius_(radius),
        scintilator_length_(scintilator_length),
        n_y_pixels_(n_y_pixels),
        n_z_pixels_(n_z_pixels),
        pixel_width_(pixel_width),
        pixel_height_(pixel_height),
        sigma_z_(sigma_z),
        sigma_dl_(sigma_dl),
        grid_center_y_(grid_center_y),
        grid_center_z_(grid_center_z),
        grid_size_y_(n_y_pixels_ * pixel_height_),
        grid_size_z_(n_z_pixels_ * pixel_width_),
        grid_ul_y_(grid_center_y_ + 0.5 * grid_size_y_),
        grid_ul_z_(grid_center_z_ - 0.5 * grid_size_z_) {
    F s_z_sq = s_z() * s_z();
    F s_dl_sq = s_dl() * s_dl();

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

  F radius() const { return radius_; }
  F scintilator_length() const { return scintilator_length_; }
  F half_scintilator_length() const { return F(0.5) * scintilator_length_; }
  F s_z() const { return sigma_z_; }
  F s_dl() const { return sigma_dl_; }

  F inv_c(int i, int j) const { return inverse_correlation_matrix_[i][j]; }

  event<F> to_projection_space_tan(const ImageSpaceEventTan<F>& ev) {
    F z_u = ev.z + (radius() - ev.y) * ev.tan;
    F z_d = ev.z - (radius() + ev.y) * ev.tan;
    F dl = -F(2.0) * ev.y * sqrt(ev.tan * ev.tan + F(1.0));
    return event<F>(z_u, z_d, dl);
  }

  event<F> to_projection_space_angle(const ImageSpaceEventAngle<F>& img_event) {
    return to_angle(to_projection_space_angle(img_event));
  }

  ImageSpaceEventTan<F> from_projection_space_tan(const event<F>& ev) {
    F t = event_tan(ev.z_u, ev.z_d, radius());
    F y = event_y(ev.dl, t);
    F z = event_z(ev.z_u, ev.z_d, y, t);
    return ImageSpaceEventTan<F>(y, z, t);
  }

  ImageSpaceEventAngle<F> from_projection_space_angle(const event<F>& ev) {
    return to_angle(from_projection_space_tan(ev));
  }

  Point pixel_center(int i, int j) {
    return std::make_pair<F>(
        grid_ul_y_ - i * pixel_height_ - 0.5 * pixel_height_,
        grid_ul_z_ + j * pixel_width_ + 0.5 * pixel_width_);
  }

  Point pixel_center(Pixel pix) { return pixel_center(pix.first, pix.second); }

  Pixel pixel_location(F y, F z) {
    return std::make_pair<int>(floor((grid_ul_y_ - y) / pixel_height_),
                               floor((z - grid_ul_z_) / pixel_width_));
  }

  Pixel pixel_location(Point p) { return pixel_location(p.first, p.second); }

  F sensitivity(F y, F z) {
    F L_plus = (half_scintilator_length() + z);
    F L_minus = (half_scintilator_length() - z);
    F R_plus = radius() + y;
    F R_minus = radius() - y;

    return INVERSE_PI *
           (std::atan(std::min(L_minus / R_minus, L_plus / R_plus)) -
            std::atan(std::max(-L_plus / R_minus, -L_minus / R_plus)));
  }

  F get_sigma_dl() { return sigma_dl_; }
  F get_sigma_z() { return sigma_z_; }

 private:
  const F radius_;
  const F scintilator_length_;
  const int n_y_pixels_;
  const int n_z_pixels_;
  const F pixel_width_;
  const F pixel_height_;
  const F sigma_z_;
  const F sigma_dl_;
  const F grid_center_y_;
  const F grid_center_z_;

  F grid_size_y_;
  F grid_size_z_;
  F grid_ul_y_;
  F grid_ul_z_;

  F inverse_correlation_matrix_[3][3];
};

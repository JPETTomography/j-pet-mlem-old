#pragma once

#include "event.h"

#include "geometry/pixel.h"
#include "geometry/point.h"
#include "util/compat.h"

/// Class responsible for the Strip detector together with the pixel grid
/// inside.
template <typename FType = double> class StripDetector {
 public:
  typedef FType F;
  typedef FType FVec[3];
  typedef ::Pixel<> Pixel;
  typedef ::Point<F> Point;

 private:
  static constexpr const F INVERSE_PI = F(M_1_PI);

 public:
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
        inverse_correlation_matrix_diag{ 1 / (sigma_z * sigma_z),
                                         1 / (sigma_z * sigma_z),
                                         1 / (sigma_dl * sigma_dl) } {}

  F half_scintilator_length() const { return F(0.5) * scintilator_length; }

  Event<F> to_projection_space_tan(const ImageSpaceEventTan<F>& ev) {
    F z_u = ev.z + (radius - ev.y) * ev.tan;
    F z_d = ev.z - (radius + ev.y) * ev.tan;
    F dl = -2 * ev.y * sqrt(ev.tan * ev.tan + 1);
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
    return Point(grid_ul_y - i * pixel_height - F(0.5) * pixel_height,
                 grid_ul_z + j * pixel_width + F(0.5) * pixel_width);
  }

  Point pixel_center(Pixel pixel) { return pixel_center(pixel.x, pixel.y); }

  Pixel pixel_location(F y, F z) {
    return Pixel(compat::floor((grid_ul_y - y) / pixel_height),
                 compat::floor((z - grid_ul_z) / pixel_width));
  }

  Pixel pixel_location(Point p) { return pixel_location(p.x, p.y); }

  F sensitivity(F y, F z) {
    F L_plus = (half_scintilator_length() + z);
    F L_minus = (half_scintilator_length() - z);
    F R_plus = radius + y;
    F R_minus = radius - y;

    return INVERSE_PI *
           (compat::atan(compat::min(L_minus / R_minus, L_plus / R_plus)) -
            compat::atan(compat::max(-L_plus / R_minus, -L_minus / R_plus)));
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
  const FVec inverse_correlation_matrix_diag;

 private:
  F grid_size_y;
  F grid_size_z;
  F grid_ul_y;
  F grid_ul_z;
};

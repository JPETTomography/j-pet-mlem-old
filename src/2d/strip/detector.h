/// \page cmd_2d_strip 2D Strip
/// \brief 2D Strip Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_2d_strip_phantom
/// - \subpage cmd_2d_strip_reconstruction

#pragma once

#if !__CUDACC__
#include <random>
#endif

#include "event.h"

#include "2d/geometry/pixel.h"
#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Strip {

/// Strip-detector made of two strip scintillators and pixel grid
template <typename FType = double> class Detector {
 public:
  using F = FType;
  using Pixel = PET2D::Pixel<>;
  using Point = PET2D::Point<F>;

  /// Creates strip-detector with given parameters.
  Detector(F radius,               //< radius of strip-detector along y-axis
           F scintillator_length,  //< lenght of strip along z-axis
           int n_y_pixels,         //< number of pixels along y-axis
           int n_z_pixels,         //< number of pixels along z-axis
           F pixel_height,         //< pixel size along y-axis
           F pixel_width,          //< pixel size along z-axis
           F sigma_z,              //< sigma z
           F sigma_dl,             //< sigma dl
           F center_y = 0,         //< center of detector y coordinate
           F center_z = 0          //< center of detector z coordinate
           )
      : radius(radius),
        scintillator_length(scintillator_length),
        n_y_pixels(n_y_pixels),
        n_z_pixels(n_z_pixels),
        total_n_pixels(n_y_pixels * n_z_pixels),
        pixel_width(pixel_width),
        pixel_height(pixel_height),
        sigma_z(sigma_z),
        sigma_dl(sigma_dl),
        center_y(center_y),
        center_z(center_z),
        size_y(n_y_pixels * pixel_height),
        size_z(n_z_pixels * pixel_width),
        tl_y(center_y + size_y / 2),
        tl_z(center_z - size_z / 2),
        tl_y_half_h(tl_y - pixel_height / 2),
        tl_z_half_w(tl_z + pixel_width / 2),
        half_scintillator_length(scintillator_length / 2) {}

#ifndef __CUDACC__
  /// Creates strip-detector determining pixel size from given detector
  /// dimensions and number of pixels.
  Detector(F radius,               //< radius of strip-detector along y-axis
           F scintillator_length,  //< lenght of strip along z-axis
           int n_y_pixels,         //< number of pixels along y-axis
           int n_z_pixels,         //< number of pixels along z-axis
           F sigma_z,              //< sigma z
           F sigma_dl              //< sigma dl
           )
      : Detector(radius,
                 scintillator_length,
                 n_y_pixels,
                 n_z_pixels,
                 2 * radius / n_y_pixels,
                 scintillator_length / n_z_pixels,
                 sigma_z,
                 sigma_dl,
                 0,
                 0) {}

  /// Copy-construct strip-detector from other detector using different float
  /// representation.
  template <typename OtherFType>
  Detector(const Detector<OtherFType>& other)
      : Detector(other.radius,
                 other.scintillator_length,
                 other.n_y_pixels,
                 other.n_z_pixels,
                 other.pixel_width,
                 other.pixel_height,
                 other.sigma_z,
                 other.sigma_dl,
                 other.center_y,
                 other.center_z) {}
#endif

  /// Convert image space event tangent to projection space.
  Event<F> to_projection_space_tan(
      const ImageSpaceEventTan<F>& is_event) const {
    F z_u = is_event.z + (radius - is_event.y) * is_event.tan;
    F z_d = is_event.z - (radius + is_event.y) * is_event.tan;
    F dl = -2 * is_event.y * sqrt(is_event.tan * is_event.tan + 1);
    return Event<F>(z_u, z_d, dl);
  }

  /// Convert image space event angle to projection space.
  Event<F> to_projection_space_angle(
      const ImageSpaceEventAngle<F>& is_ea) const {
    return to_projection_space_tan(is_ea.to_tan());
  }

  /// Convert project space event to image space event tangent.
  ImageSpaceEventTan<F> from_projection_space_tan(const Event<F>& event) const {
    F tan, y, z;
    event.transform(radius, tan, y, z);
    return ImageSpaceEventTan<F>(y, z, tan);
  }

  /// Convert project space event to image space event angle.
  ImageSpaceEventAngle<F> from_projection_space_angle(
      const Event<F>& event) const {
    return from_projection_space_tan(event).to_angle();
  }

  /// Returns pixel center point for given pixel.
  _ Point pixel_center(Pixel p) const {
    return Point(tl_z_half_w + p.x * pixel_width,
                 tl_y_half_h - p.y * pixel_height);
  }

  /// Returns pixel at given point.
  _ Pixel pixel_at(Point p) const {
    return Pixel((p.x - tl_z) / pixel_width, (tl_y - p.y) / pixel_height);
  }

  /// Returns sensitiviy of detector at given point.
  _ F sensitivity(Point p) const {
    F L_plus = half_scintillator_length + p.x;
    F L_minus = half_scintillator_length - p.x;
    F R_plus = radius + p.y;
    F R_minus = radius - p.y;

    if (R_plus <= 0) {
      R_plus = 1;
    }
    if (R_minus <= 0) {
      R_minus = 1;
    }

#if __TEST_SENSITIVITY__
    printf("L_p: %f L_m: %f\n", L_plus, L_minus);
    printf("R_p: %f R_m: %f\n", R_plus, R_minus);
    printf("FIRST: %f SECOND: %f\n",
           compat::atan(compat::min(L_minus / R_minus, L_plus / R_plus)),
           compat::atan(compat::max(-L_plus / R_minus, -L_minus / R_plus)));
#endif

    return F(M_1_PI) *
           (compat::atan(compat::min(L_minus / R_minus, L_plus / R_plus)) -
            compat::atan(compat::max(-L_plus / R_minus, -L_minus / R_plus)));
  }

  /// Returns sensitivity approximation at given pixel.
  _ F pixel_sensitivity(Pixel p) const {

    Point point = this->pixel_center(p);
    Point ur(pixel_width / 2, pixel_height / 2);
    Point ul(-pixel_width / 2, pixel_height / 2);

    return this->sensitivity(point) / 3 +       // center
           this->sensitivity(point + ur) / 6 +  // top-right
           this->sensitivity(point - ur) / 6 +  // bottom-left
           this->sensitivity(point + ul) / 6 +  // top-left
           this->sensitivity(point - ul) / 6;   // bottom-right
  }

  /// Checks if pixel is in the detector.
  bool contains_pixel(Pixel p) {
    return p.x >= 0 && p.x < (this->n_z_pixels) &&  //
           p.y >= 0 && p.y < (this->n_y_pixels);
  }

#if !__CUDACC__
  template <typename G>
  std::pair<Event<F>, bool> detect_event(const ImageSpaceEventAngle<F> is_event,
                                         G& gen) {

    Event<F> ps_event = to_projection_space_angle(is_event);

    F z_u, z_d, dl;
    z_u = ps_event.z_u;
    z_d = ps_event.z_d;
    dl = ps_event.dl;

    std::normal_distribution<F> normal_dist_dz(0, sigma_z);
    std::normal_distribution<F> normal_dist_dl(0, sigma_dl);

    z_u += normal_dist_dz(gen);
    z_d += normal_dist_dz(gen);
    dl += normal_dist_dl(gen);

    if (std::abs(z_u) < scintillator_length / 2 &&
        std::abs(z_d) < scintillator_length / 2) {

      Event<F> event(z_u, z_d, dl);
      return std::make_pair(event, true);
    } else
      return std::make_pair(Event<F>(0, 0, 0), false);
  }
#endif

  const F radius;
  const F scintillator_length;
  const int n_y_pixels;
  const int n_z_pixels;
  const int total_n_pixels;
  const F pixel_width;
  const F pixel_height;
  const F sigma_z;
  const F sigma_dl;
  const F center_y;
  const F center_z;
  const F size_y;
  const F size_z;
  const F tl_y;
  const F tl_z;
  const F tl_y_half_h;
  const F tl_z_half_w;

 private:
  const F half_scintillator_length;
};
}  // Strip
}  // PET2D

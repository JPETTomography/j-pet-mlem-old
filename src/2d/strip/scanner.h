#pragma once

#if !__CUDACC__
#include <random>
#endif

#include "event.h"

#include "2d/geometry/pixel.h"
#include "2d/geometry/point.h"
#include "util/cuda/compat.h"

/// Two-dimensional PET
namespace PET2D {
/// Two-dimensional strip PET
namespace Strip {

/// Strip-scanner made of two strip detectors (scintillators) and pixel grid

/// Represents scanner made of two scintillator strips and pixel grid
///
/// \image html detector_frame.pdf.png

template <typename FType, typename SType> class Scanner {
 public:
  using F = FType;
  using S = SType;
  using Pixel = PET2D::Pixel<S>;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Response = PET2D::Strip::Event<F>;
  using FullResponse = Response;
  using Event = PET2D::Event<F>;

  /// Creates strip-scanner with given parameters.
  Scanner(F radius,               //< radius of strip-scanner along y-axis
          F scintillator_length,  //< lenght of strip along z-axis
          int n_y_pixels,         //< number of pixels along y-axis
          int n_z_pixels,         //< number of pixels along z-axis
          F pixel_height,         //< pixel size along y-axis
          F pixel_width,          //< pixel size along z-axis
          F sigma_z,              //< sigma z
          F sigma_dl,             //< sigma dl
          F center_y = 0,         //< center of scanner y coordinate
          F center_z = 0          //< center of scanner z coordinate
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
  /// Creates strip-scanner determining pixel size from given scanner
  /// dimensions and number of pixels.
  Scanner(F radius,               //< radius of strip-scanner along y-axis
          F scintillator_length,  //< lenght of strip along z-axis
          int n_y_pixels,         //< number of pixels along y-axis
          int n_z_pixels,         //< number of pixels along z-axis
          F sigma_z,              //< sigma z
          F sigma_dl              //< sigma dl
          )
      : Scanner(radius,
                scintillator_length,
                n_y_pixels,
                n_z_pixels,
                2 * radius / n_y_pixels,
                scintillator_length / n_z_pixels,
                sigma_z,
                sigma_dl,
                0,
                0) {}

  /// Copy-construct strip-scanner from other scanner using different float
  /// representation.
  template <typename OtherFType, typename OtherSType>
  Scanner(const Scanner<OtherFType, OtherSType>& other)
      : Scanner(other.radius,
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
  Response to_projection_space_tan(
      const ImageSpaceEventTan<F>& is_event) const {
    F z_u = is_event.z + (radius - is_event.y) * is_event.tan;
    F z_d = is_event.z - (radius + is_event.y) * is_event.tan;
    F dl = -2 * is_event.y * sqrt(is_event.tan * is_event.tan + 1);
    return Response(z_u, z_d, dl);
  }

  /// Convert image space event angle to projection space.
  Response to_projection_space_angle(
      const ImageSpaceEventAngle<F>& is_ea) const {
    return to_projection_space_tan(is_ea.to_tan());
  }

  /// Convert project space event to image space event tangent.
  ImageSpaceEventTan<F> from_projection_space_tan(const Response& event) const {
    F tan, y, z;
    event.transform(radius, tan, y, z);
    return ImageSpaceEventTan<F>(y, z, tan);
  }

  /// Convert project space event to image space event angle.
  ImageSpaceEventAngle<F> from_projection_space_angle(
      const Response& event) const {
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

  /// Returns sensitiviy of scanner at given point.
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
    Vector ur(pixel_width / 2, pixel_height / 2);
    Vector ul(-pixel_width / 2, pixel_height / 2);

    return this->sensitivity(point) / 3 +       // center
           this->sensitivity(point + ur) / 6 +  // top-right
           this->sensitivity(point - ur) / 6 +  // bottom-left
           this->sensitivity(point + ul) / 6 +  // top-left
           this->sensitivity(point - ul) / 6;   // bottom-right
  }

  /// Checks if pixel is in the scanner.
  bool contains_pixel(Pixel p) {
    return p.x >= 0 && p.x < (this->n_z_pixels) &&  //
           p.y >= 0 && p.y < (this->n_y_pixels);
  }

#if !__CUDACC__

  template <class RNG>
  _ short exact_detect(
      RNG& rng,               ///< random number generator
      const Event& event,     ///< event to be detected
      FullResponse& response  ///< scanner response (LOR+length)
      ) const {

    FullResponse ps_event = to_projection_space_angle(event);

    F z_u, z_d, dl;
    z_u = ps_event.z_u;
    z_d = ps_event.z_d;
    dl = ps_event.dl;

    std::normal_distribution<F> normal_dist_dz(0, sigma_z);
    std::normal_distribution<F> normal_dist_dl(0, sigma_dl);

    z_u += normal_dist_dz(rng);
    z_d += normal_dist_dz(rng);
    dl += normal_dist_dl(rng);

    if (std::abs(z_u) < scintillator_length / 2 &&
        std::abs(z_d) < scintillator_length / 2) {
      response = Response(z_u, z_d, dl);
      return 2;
    } else
      return 0;
  }

  /// This function must be present to make the compatible interface with
  /// PhantomMonteCarlo
  template <class RNG, class AcceptanceModel>
  _ short exact_detect(RNG& rng,                ///< random number generator
                       AcceptanceModel& model,  ///< (unused)
                       const Event& event,      ///< event to be detected
                       Response& response  ///< scanner response (LOR+length)
                       ) const {
    (void)model;  // mark as unused
    return exact_detect(rng, event, response);
  }

  template <class RNG>
  _ short detect(RNG& rng,            ///< random number generator
                 const Event& event,  ///< event to be detected
                 Response& response   ///< scanner response (LOR+length)
                 ) const {
    return exact_detect(rng, event, response);
  }

  template <class RNG>
  std::pair<Response, bool> detect_event(const Event is_event, RNG& rng) {
    Response response;
    if (detect(rng, is_event, response) > 1)
      return std::make_pair(response, true);
    else
      return std::make_pair(response, false);
  }

  Response response_wo_error(const FullResponse& full_response) {
    return full_response;
  }

  template <class RNG>
  Response response_w_error(RNG& rng, const FullResponse& full_response) {
    std::normal_distribution<F> dist_z(0, sigma_z);
    std::normal_distribution<F> dist_dl(0, sigma_dl);
    Response response = response_wo_error(full_response);

    response.z_u += dist_z(rng);
    response.z_d += dist_z(rng);
    response.dl += dist_dl(rng);

    return response;
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

#pragma once

#include "cmdline.h"
#include "2d_strip/strip_detector.h"

void add_detector_options(cmdline::parser& parser);
void add_reconstruction_options(cmdline::parser& parser);
void add_phantom_options(cmdline::parser& parser);

template <typename F>
StripDetector<F> strip_detector_from_options(const cmdline::parser& parser) {

  auto R_distance = parser.get<double>("r-distance");
  auto scintillator_length = parser.get<double>("s-length");
  auto pixel_size = parser.get<double>("p-size");

  auto sigma_z = parser.get<double>("s-z");
  auto sigma_dl = parser.get<double>("s-dl");

  int n_z_pixels;
  int n_y_pixels;
  if (parser.exist("n-pixels")) {
    n_z_pixels = parser.get<int>("n-pixels");
    n_y_pixels = parser.get<int>("n-pixels");
  } else {
    if (parser.exist("n-z-pixels"))
      n_z_pixels = parser.get<int>("n-z-pixels");
    else
      n_z_pixels = (int)std::ceil(scintillator_length / pixel_size);

    if (parser.exist("n-y-pixels"))
      n_y_pixels = parser.get<int>("nyz-pixels");
    else
      n_y_pixels = (int)std::ceil(2 * R_distance / pixel_size);
  }

  return StripDetector<F>(R_distance,
                          scintillator_length,
                          n_y_pixels,
                          n_z_pixels,
                          pixel_size,
                          pixel_size,
                          sigma_z,
                          sigma_dl);
}

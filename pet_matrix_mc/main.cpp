// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detectorilators

#include <random>
#include <iostream>
#include <cmdline.h>

#include "detector_ring.h"

template <typename F> F deg(F rad) { return rad * 180/M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI/180; }

namespace cmdline {
  namespace detail {
    template <> inline std::string readable_typename<size_t>() { return "size"; }
    template <> inline std::string readable_typename<double>() { return "float"; }

    template <>
    inline std::string default_value<double>(double def) {
      if (def == 0.) return "auto";
      return detail::lexical_cast<std::string>(def);
    }
  }
}

template <typename F = double>
class always_accept {
public:
  always_accept() {}
  bool operator () (F) { return true; }
};

int main(int argc, char *argv[]) {

  cmdline::parser cl;

  cl.add<size_t>("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>("n-emissions", 'e', "emissions per pixel",               false, 1);
  cl.add<double>("radious",     'r', "inner detector ring radious",       false);
  cl.add<double>("s-pixel",     'p', "pixel size",                        false);
  cl.add<double>("w-detector",  'w', "detector width",                    false);
  cl.add<double>("h-detector",  'h', "detector height",                   false);

  cl.parse_check(argc, argv);

  auto n_pixels    = cl.get<size_t>("n-pixels");
  auto n_detectors = cl.get<size_t>("n-detectors");
  auto n_emissions = cl.get<size_t>("n-emissions");
  auto radious     = cl.get<double>("radious");
  auto s_pixel     = cl.get<double>("s-pixel");
  auto w_detector  = cl.get<double>("w-detector");
  auto h_detector  = cl.get<double>("h-detector");

  // automatic pixel size
  if (radious == 0.) {
    if (cl.get<double>("s-pixel") == 0.) {
      radious = sqrt(2.);
    } else {
      radious = sqrt(s_pixel * n_pixels);
    }
    std::cerr << "--radious=" << radious << std::endl;
  }

  // automatic radious
  if (s_pixel == 0.) {
    if (cl.get<double>("radious") == 0.) {
      s_pixel = 2./n_pixels;
    } else {
      s_pixel = radious*radious / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (w_detector == 0.) {
    w_detector = 2 * M_PI * .9 * radious / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (h_detector == 0.) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  detector_ring<> dr(n_detectors, n_pixels, s_pixel, radious, w_detector, h_detector);
  dr.matrix_mc(gen, always_accept<>(), n_emissions, true, true);
  
  auto pixel_max = 0;
  auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
  for (auto y = 0; y < n_pixels; ++y) {
    for (auto x = 0; x < n_pixels; ++x) {
      auto hits = dr.hits(x, y);
      pixel_max = std::max(pixel_max, hits);
      pixel_min = std::min(pixel_min, hits);
    }
  }

  std::cerr << "Done." << std::endl;
  std::cerr
    << "Non zero LORs: "
    << dr.non_zero_lors()
    << '/'
    << dr.lors()
    << std::endl;
  std::cerr
    << "Min hits: "
    << pixel_min
    << std::endl;
  std::cerr
    << "Max hits: "
    << pixel_max
    << std::endl;
  std::cerr << "Press Enter." << std::endl;
  while(getc(stdin) != '\n');

  std::cout << "P2" << std::endl;
  std::cout << "# pet_matrix.pgm" << std::endl;
  std::cout << n_pixels << ' ' << n_pixels << std::endl;
  std::cout << pixel_max << std::endl;
  for (auto y = 0; y < n_pixels; ++y) {
    for (auto x = 0; x < n_pixels; ++x) {
      std::cout << dr.hits(x, y) << ' ';
    }
    std::cout << std::endl;
  }

  return 0;
}

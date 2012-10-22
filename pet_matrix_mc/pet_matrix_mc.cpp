// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detectorilators

#include <random>
#include <iostream>
#include <cmdline.h>

#include "circle.h"
#include "event.h"
#include "detector.h"
#include "reconstruction.h"

template <typename F> F deg(F rad) { return rad * 180/M_PI; }
template <typename F> F rad(F deg) { return deg * M_PI/180; }

typedef std::map<std::pair<int, int>, std::vector<int>> lor_map;

template <>
inline std::string cmdline::detail::readable_typename<size_t>() { return "size"; }
template <>
inline std::string cmdline::detail::readable_typename<double>() { return "float"; }

template <>
inline std::string cmdline::detail::default_value<double>(double def) {
  if (def == 0.) return "auto";
  return detail::lexical_cast<std::string>(def);
}

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

#if 1
  lor_map mc;
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> one_dis(0, 1);
  std::uniform_real_distribution<> phi_dis(0, M_PI);
  circle<> c_inner(radious);
  circle<> c_outer(radious+h_detector);

  // iterating only triangular matrix,
  // being upper right part or whole system matrix
  for (auto y = 0; y < n_pixels/2; ++y)
    for (auto x = 0; x <= y; ++x)
      for (auto n = 0; n < n_emissions; ++n) {
        auto rx = x + one_dis(gen);
        auto ry = y + one_dis(gen);
        // ensure we are within a triangle
        if (rx > ry) continue;
        // random point within a pixel
        decltype(c_inner)::event_type e(rx, ry, phi_dis(gen));
        // secant for p and phi
        auto s_inner = c_inner.secant(e);
        auto s_outer = c_outer.secant(e);

        std::cout << '(' << e.x << ',' << e.y << ')'
                  << ' ' << deg(e.phi)
                  << " s1 = (" << s_inner.first.x  << ',' << s_inner.first.y  << ')'
                  << ' ' << deg(atan2(s_inner.first.y, s_inner.first.x))
                  << " s2 = (" << s_inner.second.x << ',' << s_inner.second.y << ')'
                  << ' ' << deg(atan2(s_inner.second.y, s_inner.second.x))
                  << std::endl;
      }
#else
  Ring2DDetector<int, double> detector(n_detectors, n_pixels);
  ProbabilityMatrix<int, double> probability_matrix(detector);

  const int max_lors = detector.n_detectors()*detector.n_detectors();

  double *p = new double[max_lors];

  for (int ipix = 0; ipix < probability_matrix.octant_size(); ++ipix) {
    auto pix = probability_matrix.octant(ipix);

    for (int i = 0; i<max_lors; ++i) p[i]=0.0;

    fill_pixel(pix, detector.n_detectors(), p, n_emissions);

    auto pixel_row = row_from_array(pix, p,detector.n_detectors(), max_lors);

    probability_matrix.push_back_row(pixel_row);
    std::cout << *pixel_row;
  }

  FILE *fout = fopen("prob.bin", "w");
  probability_matrix.fwrite(fout);
  fclose(fout);
#endif

  return 0;
}

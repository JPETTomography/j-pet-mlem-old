// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detector scintilators.

#include <random>
#include <iostream>

#include <cmdline.h>
#include "cmdline_types.h"

#include "detector_ring.h"
#include "model.h"
#include "png_writer.h"
#include "svg_ostream.h"

#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[]) {

try {
  cmdline::parser cl;
  cl.footer("matrix_file ...");

#if _OPENMP
  cl.add<size_t>     ("n-threads",   't', "number of OpenMP threads",          false);
#endif
  cl.add<size_t>     ("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>     ("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>     ("n-emissions", 'e', "emissions per pixel",               false, 1);
  cl.add<double>     ("radious",     'r', "inner detector ring radious",       false);
  cl.add<double>     ("s-pixel",     'p', "pixel size",                        false);
  cl.add<double>     ("w-detector",  'w', "detector width",                    false);
  cl.add<double>     ("h-detector",  'h', "detector height",                   false);
  cl.add<std::string>("model",       'm', "acceptance model",                  false,
                      "scintilator", cmdline::oneof<std::string>("always", "scintilator"));
  cl.add<double>     ("acceptance",  'a', "acceptance probability factor",     false, 10.);
  cl.add             ("stats",         0, "show stats");
  cl.add             ("wait",          0, "wait before exit");
  cl.add<ssize_t>    ("lor",         'l', "select lor to output to a file",    false, -1);
  cl.add<std::string>("output",      'o', "output binary triangular sparse system matrix", false);
  cl.add<std::string>("png",           0, "output PNG with hit/system matrix", false);
  cl.add<std::string>("svg",           0, "output SVG detector ring geometry", false);
  cl.add<std::mt19937::result_type>
                     ("seed",        's', "random number generator seed",      false);

  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif

  auto n_pixels    = cl.get<size_t>("n-pixels");
  auto n_detectors = cl.get<size_t>("n-detectors");
  auto n_emissions = cl.get<size_t>("n-emissions");
  auto radious     = cl.get<double>("radious");
  auto s_pixel     = cl.get<double>("s-pixel");
  auto w_detector  = cl.get<double>("w-detector");
  auto h_detector  = cl.get<double>("h-detector");

  // automatic pixel size
  if (!cl.exist("radious")) {
    if (!cl.exist("s-pixel")) {
      radious = M_SQRT2; // exact result
    } else {
      radious = s_pixel * n_pixels / M_SQRT2;
    }
    std::cerr << "--radious=" << radious << std::endl;
  }

  // automatic radious
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radious")) {
      s_pixel = 2. / n_pixels; // exact result
    } else {
      s_pixel = M_SQRT2 * radious / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radious / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<std::mt19937::result_type>("seed"));
  }

  detector_ring<> dr(n_detectors, n_pixels, s_pixel, radious, w_detector, h_detector);

  for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
    ibstream in(*fn, std::ios::binary);
    in >> dr;
  }

  if (cl.get<std::string>("model") == "always")
    dr.matrix_mc(gen, always_accept<>(), n_emissions);
  if (cl.get<std::string>("model") == "scintilator")
    dr.matrix_mc(
      gen,
      scintilator_accept<>(cl.get<double>("acceptance")),
      n_emissions);

  auto pixel_max = 0;
  auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
  auto lor       = cl.get<ssize_t>("lor");

  if (cl.exist("stats") || cl.exist("png"))
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        auto hits = lor > 0 ? dr.matrix(lor, x, y) : dr.hits(x, y);
        pixel_max = std::max(pixel_max, hits);
        pixel_min = std::min(pixel_min, hits);
      }
    }

  // show stats if requested
  if (cl.exist("stats")) {
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
  }

  if (cl.exist("output")) {
    obstream out(cl.get<std::string>("output"), std::ios::binary | std::ios::trunc);
    out << dr;
  }

  if (cl.exist("png")) {
    png_writer png(cl.get<std::string>("png"));
    dr.output_bitmap(png, pixel_max, lor);
  }

  if (cl.exist("svg")) {
    svg_ostream<> svg(cl.get<std::string>("svg"),
      radious+h_detector, radious+h_detector,
      1024., 1024.);
    svg << dr;
    if (cl.exist("png")) {
      svg.link_image(cl.get<std::string>("png"),
        -(s_pixel*n_pixels)/2., -(s_pixel*n_pixels)/2.,
          s_pixel*n_pixels, s_pixel*n_pixels);
    }
  }

  if (cl.exist("wait")) {
    std::cerr << "Press Enter." << std::endl;
    while(getc(stdin) != '\n');
  }

  return 0;

} catch(std::string &ex) {
  std::cerr << "error: " << ex << std::endl;
} catch(const char *ex) {
  std::cerr << "error: " << ex << std::endl;
}

}

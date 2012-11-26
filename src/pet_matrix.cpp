// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detector scintilators.

#include <iostream>
#include <random>

#include <cmdline.h>
#include "cmdline_types.h"

#include "random.h"
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
  cl.add<size_t>     ("n-threads",   't', "number of OpenMP threads",          false, 0, false);
#endif
  cl.add<size_t>     ("n-pixels",    'n', "number of pixels in one dimension", false, 256);
  cl.add<size_t>     ("n-detectors", 'd', "number of ring detectors",          false, 64);
  cl.add<size_t>     ("n-emissions", 'e', "emissions per pixel",               false, 0);
  cl.add<double>     ("radious",     'r', "inner detector ring radious",       false, 0);
  cl.add<double>     ("s-pixel",     'p', "pixel size",                        false);
  cl.add<double>     ("w-detector",  'w', "detector width",                    false);
  cl.add<double>     ("h-detector",  'h', "detector height",                   false);
  cl.add<std::string>("model",       'm', "acceptance model",                  false,
                      "scintilator", cmdline::oneof<std::string>("always", "scintilator"));
  cl.add<double>     ("acceptance",  'a', "acceptance probability factor",     false, 10.);
  cl.add<tausworthe::seed_type>
                     ("seed",        's', "random number generator seed",      false);
  cl.add<std::string>("output",      'o', "output binary triangular/full sparse system matrix", false, "", false);
  cl.add             ("full",        'f', "output full non-triangular sparse system matrix");

  // visual debugging params
  cl.add<std::string>("png",           0, "output lor to png",            false, "", false);
  cl.add<ssize_t>    ("from",          0, "lor start detector to output", false, -1, false);
  cl.add<ssize_t>    ("to",            0, "lor end detector to output",   false, -1, false);

  // printing & stats params
  cl.add             ("print",         0, "print triangular sparse system matrix");
  cl.add             ("stats",         0, "show stats");
  cl.add             ("wait",          0, "wait before exit");

  cl.parse_check(argc, argv);

  auto &n_pixels    = cl.get<size_t>("n-pixels");
  auto &n_detectors = cl.get<size_t>("n-detectors");
  auto &n_emissions = cl.get<size_t>("n-emissions");
  auto &radious     = cl.get<double>("radious");
  auto &s_pixel     = cl.get<double>("s-pixel");
  auto &w_detector  = cl.get<double>("w-detector");
  auto &h_detector  = cl.get<double>("h-detector");

  // check options
  if (cl.exist("png") && !cl.exist("from")) {
    throw("need to specify output --png option when --from is specified");
  }
  if (!cl.exist("png") && cl.exist("from")) {
    throw("need to specify --from lor when output --png option is specified");
  }

  // load config files
  for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
    std::ifstream in(*fn+".cfg");
    if (!in.is_open()) continue;
    // load except n-emissions
    auto n_prev_emissions = n_emissions;
    in >> cl;
    n_emissions = n_prev_emissions;
    break; // only one config file allowed!
  }

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<size_t>("n-threads"));
  }
#endif

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
  tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<tausworthe::seed_type>("seed"));
  }

  detector_ring<> dr(n_detectors, n_pixels, s_pixel, radious, w_detector, h_detector);

  for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
    ibstream in(*fn, std::ios::binary);
    if (!in.is_open()) throw("cannot open input file: " + *fn);
    in >> dr;
  }

  if (cl.get<std::string>("model") == "always")
    dr.matrix_mc(gen, always_accept<>(), n_emissions);
  if (cl.get<std::string>("model") == "scintilator")
    dr.matrix_mc(
      gen,
      scintilator_accept<>(cl.get<double>("acceptance")),
      n_emissions);

  // generate output
  if (cl.exist("output")) {
    auto fn = cl.get<std::string>("output");
    obstream out(fn, std::ios::binary | std::ios::trunc);
    dr.output_triangular = !cl.exist("full");
    out << dr;

    std::ofstream os(fn+".cfg", std::ios::trunc);
    os << cl;

    png_writer png(fn+".png");
    dr.output_bitmap(png);

    svg_ostream<> svg(fn+".svg",
      radious+h_detector, radious+h_detector,
      1024., 1024.);
    svg << dr;

    auto fn_sep     = fn.find_last_of("\\/");
    auto fn_wo_path = fn.substr(fn_sep != std::string::npos ? fn_sep+1 : 0);
    svg.link_image(fn_wo_path+".png",
      -(s_pixel*n_pixels)/2., -(s_pixel*n_pixels)/2.,
        s_pixel*n_pixels, s_pixel*n_pixels);
  }

  // visual debugging output
  if (cl.exist("png")) {
    detector_ring<>::lor_type lor(0, 0);
    lor.first  = cl.get<ssize_t>("from");
    if (cl.exist("to")) {
      lor.second = cl.get<ssize_t>("to");
    } else {
      lor.second = (lor.first + n_detectors / 2) % n_detectors;
    }
    png_writer png(cl.get<std::string>("png"));
    dr.output_bitmap(png, lor);
  }

  // show stats if requested
  if (cl.exist("stats")) {
    auto pixel_max = 0;
    auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        auto hits = dr.hits(x, y);
        pixel_min = std::min(pixel_min, hits);
        pixel_max = std::max(pixel_max, hits);
      }
    }
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

  if (cl.exist("print")) {
    std::cout << dr;
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

// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detector scintilators.

#define LOR_MAJOR 0

#include <iostream>
#include <random>

#include "cmdline.h"
#include "cmdline_types.h"

#include "random.h"
#include "detector_ring.h"
#if LOR_MAJOR
#include "matrix_lor_major.h"
#else
#include "matrix_pixel_major.h"
#endif
#include "pixel.h"
#include "lor.h"
#include "model.h"
#include "png_writer.h"
#include "svg_ostream.h"

#include "monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    cl.footer("matrix_file ...");

    cl.add<cmdline::string>(
        "config", 'c', "load config file", false, cmdline::string(), false);
#if _OPENMP
    cl.add<int>("n-threads", 't', "number of OpenMP threads", false, 0, false);
#endif
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 256);
    cl.add<int>("n-detectors", 'd', "number of ring detectors", false, 64);
    cl.add<int>("n-emissions", 'e', "emissions per pixel", false, 0);
    cl.add<double>("radius", 'r', "inner detector ring radius", false, 0);
    cl.add<double>("s-pixel", 'p', "pixel size", false);
    cl.add<double>("w-detector", 'w', "detector width", false);
    cl.add<double>("h-detector", 'h', "detector height", false);
    cl.add<std::string>("model",
                        'm',
                        "acceptance model",
                        false,
                        "scintilator",
                        cmdline::oneof<std::string>("always", "scintilator"));
    cl.add<double>(
        "acceptance", 'a', "acceptance probability factor", false, 10.);
    cl.add<tausworthe::seed_type>(
        "seed", 's', "random number generator seed", false, 0, false);
    cl.add<cmdline::string>(
        "output",
        'o',
        "output binary triangular/full sparse system matrix",
        false,
        cmdline::string(),
        false);
    cl.add("full", 'f', "output full non-triangular sparse system matrix");

    // visual debugging params
    cl.add<cmdline::string>(
        "png", 0, "output lor to png", false, cmdline::string(), false);
    cl.add<int>("from", 0, "lor start detector to output", false, -1, false);
    cl.add<int>("to", 0, "lor end detector to output", false, -1, false);

    // printing & stats params
    cl.add("print", 0, "print triangular sparse system matrix");
    cl.add("stats", 0, "show stats");
    cl.add("wait", 0, "wait before exit");

    cl.parse_check(argc, argv);

    auto& n_pixels = cl.get<int>("n-pixels");
    auto& n_detectors = cl.get<int>("n-detectors");
    auto& n_emissions = cl.get<int>("n-emissions");
    auto& radius = cl.get<double>("radius");
    auto& s_pixel = cl.get<double>("s-pixel");
    auto& w_detector = cl.get<double>("w-detector");
    auto& h_detector = cl.get<double>("h-detector");

    // check options
    if (cl.exist("png") && !cl.exist("from")) {
      throw("need to specify output --png option when --from is specified");
    }
    if (!cl.exist("png") && cl.exist("from")) {
      throw("need to specify --from lor when output --png option is specified");
    }

    // load config file
    if (cl.exist("config")) {
      std::ifstream in(cl.get<cmdline::string>("config"));
      if (!in.is_open()) {
        throw("cannot open input config file: " +
              cl.get<cmdline::string>("config"));
      }
      // load except n-emissions
      auto n_prev_emissions = n_emissions;
      in >> cl;
      n_emissions = n_prev_emissions;
    } else
      // load config files accompanying matrix files
      for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
        auto fn_sep = fn->find_last_of("\\/");
        auto fn_ext = fn->find_last_of(".");
        auto fn_wo_ext =
            fn->substr(0,
                       fn_ext != std::string::npos &&
                       (fn_sep == std::string::npos || fn_sep < fn_ext)
                           ? fn_ext : std::string::npos);
        std::ifstream in(fn_wo_ext + ".cfg");
        if (!in.is_open())
          continue;
        // load except n-emissions
        auto n_prev_emissions = n_emissions;
        in >> cl;
        n_emissions = n_prev_emissions;
        break;  // only one config file allowed!
      }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    // automatic pixel size
    if (!cl.exist("radius")) {
      if (!cl.exist("s-pixel")) {
        radius = M_SQRT2;  // exact result
      } else {
        radius = s_pixel * n_pixels / M_SQRT2;
      }
      std::cerr << "--radius=" << radius << std::endl;
    }

    // automatic radius
    if (!cl.exist("s-pixel")) {
      if (!cl.exist("radius")) {
        s_pixel = 2. / n_pixels;  // exact result
      } else {
        s_pixel = M_SQRT2 * radius / n_pixels;
      }
      std::cerr << "--s-pixel=" << s_pixel << std::endl;
    }

    // automatic detector size
    if (!cl.exist("w-detector")) {
      w_detector = 2 * M_PI * .9 * radius / n_detectors;
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

    DetectorRing<> dr(n_detectors, radius, w_detector, h_detector);
#if LOR_MAJOR
    MatrixLORMajor<Pixel<>, LOR<>> matrix(n_pixels, n_detectors);
#else
    MatrixPixelMajor<Pixel<>, LOR<>> matrix(n_pixels, n_detectors);
#endif

    MonteCarlo<decltype(dr), decltype(matrix)> monte_carlo(dr, matrix, s_pixel);

    for (auto fn = cl.rest().begin(); fn != cl.rest().end(); ++fn) {
      ibstream in(*fn, std::ios::binary);
      if (!in.is_open())
        throw("cannot open input file: " + *fn);
      try {
        decltype(matrix) ::SparseMatrix sparse_matrix(in);
        sparse_matrix.sort_by_pixel();
        matrix << sparse_matrix;
      }
      catch (std::string & ex) {
        throw(ex + ": " + *fn);
      }
      catch (const char * ex) {
        throw(std::string(ex) + ": " + *fn);
      }
    }

    if (cl.get<std::string>("model") == "always")
      monte_carlo(gen, AlwaysAccept<>(), n_emissions);
    if (cl.get<std::string>("model") == "scintilator")
      monte_carlo(
          gen, ScintilatorAccept<>(cl.get<double>("acceptance")), n_emissions);

    auto sparse_matrix = matrix.to_sparse();

    // generate output
    if (cl.exist("output")) {
      auto fn = cl.get<cmdline::string>("output");
      auto fn_sep = fn.find_last_of("\\/");
      auto fn_ext = fn.find_last_of(".");
      auto fn_wo_ext =
          fn.substr(0,
                    fn_ext != std::string::npos &&
                    (fn_sep == std::string::npos || fn_sep < fn_ext)
                        ? fn_ext : std::string::npos);
      auto fn_wo_path =
          fn_wo_ext.substr(fn_sep != std::string::npos ? fn_sep + 1 : 0);

      obstream out(fn, std::ios::binary | std::ios::trunc);
      out << sparse_matrix;

      std::ofstream os(fn_wo_ext + ".cfg", std::ios::trunc);
      os << cl;

      png_writer png(fn_wo_ext + ".png");
      matrix.output_bitmap(png);

      svg_ostream<> svg(fn_wo_ext + ".svg",
                        radius + h_detector,
                        radius + h_detector,
                        1024.,
                        1024.);
      svg << dr;

      svg.link_image(fn_wo_path + ".png",
                     -(s_pixel * n_pixels) / 2.,
                     -(s_pixel * n_pixels) / 2.,
                     s_pixel * n_pixels,
                     s_pixel * n_pixels);
    }

    // visual debugging output
    if (cl.exist("png")) {
      LOR<> lor(0, 0);
      lor.first = cl.get<int>("from");
      if (cl.exist("to")) {
        lor.second = cl.get<int>("to");
      } else {
        lor.second = (lor.first + n_detectors / 2) % n_detectors;
      }
      png_writer png(cl.get<cmdline::string>("png"));
      matrix.to_sparse().output_lor_bitmap(png, lor);
    }

    // show stats if requested
    if (cl.exist("stats")) {
      auto pixel_max = 0;
      auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
      for (auto y = 0; y < n_pixels; ++y) {
        for (auto x = 0; x < n_pixels; ++x) {
          auto hits = matrix[decltype(matrix) ::Pixel(x, y)];
          pixel_min = std::min(pixel_min, hits);
          pixel_max = std::max(pixel_max, hits);
        }
      }
      std::cerr << "Non zero LORs: " << matrix.non_zero_lors() << '/'
                << dr.lors() << std::endl;
      std::cerr << "Min hits: " << pixel_min << std::endl;
      std::cerr << "Max hits: " << pixel_max << std::endl;
    }

    if (cl.exist("print")) {
      std::cout << sparse_matrix;
    }

    if (cl.exist("wait")) {
      std::cerr << "Press Enter." << std::endl;
      while (getc(stdin) != '\n')
        ;
    }

    return 0;

  }
  catch (std::string & ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char * ex) {
    std::cerr << "error: " << ex << std::endl;
  }
}

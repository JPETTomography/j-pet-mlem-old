#include "options.h"

#include <cmath>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/variant.h"

#include "util/random.h"  // util::random::tausworthe::seed_type
#include <random>         // std::mt19937::result_type

namespace PET2D {
namespace Barrel {

void add_detector_options(cmdline::parser& cl) {
  cl.add<double>("radius", 'r', "inner detector ring radius", false);
  cl.add<double>("radius2", 0, " ... 2nd ring", false);
  cl.add<double>("radius3", 0, " ... 3rd ring", false);
  cl.add<double>("radius4", 0, " ... 4th ring", false);
  cl.add<double>("rotation", 0, "ring rotation (0-1)", false);
  cl.add<double>("rotation2", 0, " ... 2nd ring", false);
  cl.add<double>("rotation3", 0, " ... 3rd ring", false);
  cl.add<double>("rotation4", 0, " ... 4th ring", false);
  cl.add<int>("n-detectors", 'd', "number of detectors in ring", false);
  cl.add<int>("n-detectors2", 0, " ... 2nd ring", false);
  cl.add<int>("n-detectors3", 0, " ... 3rd ring", false);
  cl.add<int>("n-detectors4", 0, " ... 4th ring", false);
}

void add_matrix_options(cmdline::parser& cl) {
  std::ostringstream msg;
  msg << "matrix_file ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in meters.";
  cl.footer(msg.str());

  cl.add<cmdline::path>("config",
                        'c',
                        "load config file",
                        cmdline::dontsave,
                        cmdline::path(),
                        cmdline::default_reader<cmdline::path>(),
                        cmdline::load);
  cl.add<int>("n-pixels", 'n', "number of pixels in one dimension", false, 256);
  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
  add_detector_options(cl);
  cl.add<int>("n-emissions",
              'e',
              "emissions per pixel",
              false,
              0,
              cmdline::not_from_file);
  cl.add<double>("s-pixel", 'p', "pixel size", false);
  cl.add<double>(
      "tof-step", 't', "TOF quantisation step for distance delta", false);
  cl.add<std::string>(
      "shape",
      's',
      "detector (scintillator) shape (square, circle, triangle, hexagon)",
      false,
      "square",
      cmdline::oneof<std::string>("square", "circle", "triangle", "hexagon"));
  cl.add<double>("w-detector", 'w', "detector width", false);
  cl.add<double>("h-detector", 'h', "detector height", false);
  cl.add<double>("d-detector",
                 0,
                 "inscribe detector shape into circle of given diameter",
                 false);
  cl.add<std::string>(
      "model",
      'm',
      "acceptance model (always, scintillator)",
      false,
      "scintillator",
      cmdline::oneof<std::string>("always",
                                  "scintillator",
                                  /* obsolete */ "scintilator"));
  // NOTE: this options is obsolete (use base-length instead)
  cl.add<double>("acceptance",
                 'a',
                 "acceptance probability factor",
                 cmdline::dontsave | cmdline::hidden,
                 10.);
  cl.add<double>("base-length",
                 'l',
                 "scintillator emission base length P(l)=1-e^(-1)",
                 false,
                 0.1);
  cl.add<cmdline::path>("output",
                        'o',
                        "output binary triangular/full sparse system matrix",
                        cmdline::dontsave);
  cl.add("full", 'f', "output full non-triangular sparse system matrix");

  // visual debugging params
  cl.add<cmdline::path>("png", 0, "output lor to png", cmdline::dontsave);
  cl.add<int>("from", 0, "lor start detector to output", cmdline::dontsave, -1);
  cl.add<int>("to", 0, "lor end detector to output", cmdline::dontsave, -1);
  cl.add<int>("pos", 0, "position to output", cmdline::dontsave, -1);

  // printing & stats params
  cl.add("print", 0, "print triangular sparse system matrix");
  cl.add("stats", 0, "show stats");
  cl.add("wait", 0, "wait before exit");
  cl.add("verbose", 'v', "prints the iterations information on std::out");
  cl.add<util::random::tausworthe::seed_type>(
      "seed", 'S', "random number generator seed", cmdline::dontsave);
#if HAVE_CUDA
  cl.add("gpu", 'G', "run on GPU (via CUDA)");
  cl.add<int>("cuda-blocks", 'B', "CUDA blocks", cmdline::dontsave, 64);
  cl.add<int>(
      "cuda-threads", 'W', "CUDA threads per block", cmdline::dontsave, 512);
#endif
#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#endif
}

void add_phantom_options(cmdline::parser& cl) {
  cl.footer("phantom_description");

  cl.add<cmdline::path>("config",
                        'c',
                        "load config file",
                        cmdline::dontsave,
                        cmdline::path(),
                        cmdline::default_reader<cmdline::path>(),
                        cmdline::load);
  cl.add<int>("n-pixels", 'n', "number of pixels in one dimension", false, 256);
  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
  add_detector_options(cl);
  cl.add<int>("n-emissions",
              'e',
              "emissions",
              cmdline::optional,
              0,
              cmdline::default_reader<int>(),
              cmdline::not_from_file);
  cl.add<double>("s-pixel", 'p', "pixel size", false);
  cl.add<double>(
      "tof-step", 't', "TOF quantisation step for distance delta", false);
  cl.add<std::string>(
      "shape",
      's',
      "detector (scintillator) shape (square, circle, triangle, hexagon)",
      false,
      "square",
      cmdline::oneof<std::string>("square", "circle", "triangle", "hexagon"));
  cl.add<double>("w-detector", 'w', "detector width", false);
  cl.add<double>("h-detector", 'h', "detector height", false);
  cl.add<double>("d-detector",
                 0,
                 "inscribe detector shape into circle of given diameter",
                 false);
  cl.add<std::string>(
      "model",
      'm',
      "acceptance model",
      false,
      "scintillator",
      cmdline::oneof<std::string>("always",
                                  "scintillator",
                                  /* obsolete */ "scintilator"));
  // NOTE: this options is obsolete (use base-length instead)
  cl.add<double>("acceptance",
                 'a',
                 "acceptance probability factor",
                 cmdline::dontsave | cmdline::hidden,
                 10.);
  cl.add<double>("base-length",
                 'l',
                 "scintillator emission base length P(l)=1-e^(-1)",
                 false,
                 0.1);
  cl.add<cmdline::path>(
      "output", 'o', "output lor hits for supplied phantom", cmdline::dontsave);
  cl.add("detected", 0, "collects detected emissions");

  // printing & stats params
  cl.add("verbose", 'v', "prints the iterations information on std::out");
  cl.add<std::mt19937::result_type>(
      "seed", 'S', "random number generator seed", cmdline::dontsave);
#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#endif
}

void add_reconstruction_options(cmdline::parser& cl) {
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);
  cl.add<cmdline::path>("system", 's', "system matrix file", true);
  cl.add<cmdline::path>("mean", 'm', "mean file", true);
  cl.add<cmdline::path>("output", 'o', "output reconstruction", false);
  cl.add<double>("threshold", 0, "discretisation treshold", false, 0.0);

  // additional options
  cl.add("verbose", 'v', "prints the iterations information on std::out");

#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", false);
#endif
}

void calculate_detector_options(cmdline::parser& cl) {
  // convert obsolete acceptance option to length scale
  auto& length_scale = cl.get<double>("base-length");
  if (cl.exist("acceptance") && !cl.exist("base-length")) {
    length_scale = 1.0 / cl.get<double>("acceptance");
  }
  // FIXME: fixup for spelling mistake, present in previous versions
  auto& model_name = cl.get<std::string>("model");
  if (model_name == "scintilator") {
    model_name = "scintillator";
  }

  if (cl.exist("verbose")) {
    std::cerr << "assumed:" << std::endl;
  }

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& radius = cl.get<double>("radius");
  auto& n_detectors = cl.get<int>("n-detectors");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& w_detector = cl.get<double>("w-detector");
  auto& d_detector = cl.get<double>("d-detector");
  auto& shape = cl.get<std::string>("shape");

  // automatic radius size
  if (!cl.exist("radius")) {
    if (!cl.exist("s-pixel")) {
      radius = M_SQRT2;  // exact result
    } else {
      radius = s_pixel * n_pixels / M_SQRT2;
    }
    std::cerr << "--radius=" << radius << std::endl;
  }

  // automatic pixel size
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radius")) {
      s_pixel = 2. / n_pixels;  // exact result
    } else {
      s_pixel = M_SQRT2 * radius / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  if (!cl.exist("w-detector")) {
    if (cl.exist("n-detectors")) {
      w_detector = 2 * M_PI * .9 * radius / n_detectors;
    } else if (cl.exist("d-detector")) {
      if (shape == "circle") {
        w_detector = d_detector;
      } else {
        auto mult = 1.;
        auto sides = 0.;
        if (shape == "triangle") {
          sides = 3.;
        } else if (shape == "square") {
          sides = 4.;
        } else if (shape == "hexagon") {
          sides = 6.;
          mult = 2.;
        } else {
          throw("cannot determine detector width for given shape");
        }
        w_detector = d_detector * std::sin(M_PI / sides) * mult;
      }
    }
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }

  // automatic detector size
  // NOTE: detector height will be determined per shape
  if (!cl.exist("n-detectors")) {
    if (cl.exist("d-detector")) {
      n_detectors =
          ((int)std::floor(
               M_PI / std::atan2(d_detector, 2 * radius + d_detector / 2)) /
           4) *
          4;
    } else {
      n_detectors =
          ((int)std::floor(M_PI / std::atan2(w_detector, 2 * radius)) / 4) * 4;
    }
    if (!n_detectors) {
      throw("detector width is too big for given detector ring radius");
    }
    std::cerr << "--n-detectors=" << n_detectors << std::endl;
  }
}
}
}
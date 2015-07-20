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

void add_scanner_options(cmdline::parser& cl) {
  // detector shape and dimensions
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

  // custom barrel dimensions
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
  cl.add<double>("fov-radius", 0, "field of view radius", false);
}

static void add_scanner_model_options(cmdline::parser& cl) {
  cl.add<cmdline::path>("config",
                        'c',
                        "load config file",
                        cmdline::dontsave,
                        cmdline::path(),
                        cmdline::default_reader<cmdline::path>(),
                        cmdline::load);
  cl.add<int>("n-pixels", 'n', "number of pixels in one dimension", false, 256);
  cl.add<double>("s-pixel", 'p', "pixel size", false);
  add_scanner_options(cl);
  cl.add<double>(
      "tof-step", 't', "TOF quantisation step for distance delta", false);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
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
}

void add_matrix_options(cmdline::parser& cl) {
  std::ostringstream msg;
  msg << "matrix_file ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in meters.";
  cl.footer(msg.str());

  add_scanner_model_options(cl);
  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
  cl.add<int>("n-emissions",
              'e',
              "emissions per pixel",
              false,
              0,
              cmdline::not_from_file);
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
  cl.footer("phantom_description ...");

  cl.add<cmdline::path>("config",
                        'c',
                        "load config file",
                        cmdline::dontsave,
                        cmdline::path(),
                        cmdline::default_reader<cmdline::path>(),
                        cmdline::load);
  cl.add<int>("n-pixels", 'n', "number of pixels in one dimension", false, 256);
  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
  add_scanner_options(cl);
  cl.add<int>("n-emissions",
              'e',
              "number of emissions",
              cmdline::optional,
              0,
              cmdline::default_reader<int>(),
              cmdline::not_from_file);
  cl.add<double>("s-pixel", 'p', "pixel size", false);
  cl.add<double>(
      "tof-step", 't', "TOF quantisation step for distance delta", false);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
  cl.add("bin", 0, "ouput number of hits in each lor position");

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

static void add_iteration_options(cmdline::parser& cl) {
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);
}

void add_reconstruction_options(cmdline::parser& cl) {
  cl.footer("means ...");
  add_iteration_options(cl);
  cl.add<cmdline::path>("system", 's', "system matrix file", true);
  cl.add<cmdline::path>("output", 'o', "output reconstruction", false);
  cl.add("no-sensitivity", 0, "do not correct for sensitivity");
  cl.add<double>("png-max", 0, "maximum value mapped to 255 in PNG", false, 0);

  // additional options
  cl.add("verbose", 'v', "prints the iterations information on std::out");

#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", false);
#endif
}

void add_lm_reconstruction_options(cmdline::parser& cl) {
  cl.footer("response ...");
  add_iteration_options(cl);
  cl.add<double>("length", 0, "length of the detector", false, 0.3);
  cl.add<cmdline::path>("geometry", 0, "geometry information", true);
  cl.add<cmdline::path>("system", 0, "system matrix file", false);
  cl.add<cmdline::path>("output", 'o', "output reconstruction", false);
  cl.add("graphics", 'g', "output mathematica .m graphics file", false);
  cl.add("event", 0, "event number", false, 0);
  add_scanner_model_options(cl);

  // additional options
  cl.add("verbose", 'v', "prints the iterations information on std::out");

#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#endif
}

void calculate_scanner_options(cmdline::parser& cl, int argc) {
  // check options
  if (!cl.exist("w-detector") && !cl.exist("d-detector") &&
      !cl.exist("n-detectors")) {
    if (argc == 1) {
      std::cerr << cl.usage()
                << "requires one of following options specified:" << std::endl
                << "  --w-detector, --d-detector, --n-detectors" << std::endl;
      exit(0);
    } else {
      throw("need to specify either --w-detector, --d-detector, --n-detectors");
    }
  }
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
  auto& fov_radius = cl.get<double>("fov-radius");

  // 1. Automatic radius size

  // This sets default radius size to sqrt(2), so the image space is (-1, 1).

  if (!cl.exist("radius")) {
    if (!cl.exist("s-pixel")) {
      radius = M_SQRT2;  // exact result
    } else {
      radius = s_pixel * n_pixels / M_SQRT2;
    }
    std::cerr << "--radius=" << radius << std::endl;
  }

  // 2. Set fov-radius matching radius

  if (!cl.exist("fov-radius")) {
    fov_radius = radius / M_SQRT2;
    std::cerr << "--fov-radius=" << fov_radius << std::endl;
  }

  // 3. Automatic pixel size

  if (!cl.exist("s-pixel")) {
    s_pixel = 2 * fov_radius / n_pixels;
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // 4. Automatic detector size

  // Retermine detector dimensions so it fits best given number of detectors
  // and radius.

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

  // 5. Automatic detector count

  // Detector count will be determined per shape, so we put as many as possible
  // into given radius, keeping the symmetry.

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

}  // Barrel
}  // PET2D

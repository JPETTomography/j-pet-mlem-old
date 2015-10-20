#include "options.h"

#include <cmath>

#include "common/options.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/variant.h"

namespace PET2D {
namespace Strip {

void add_scanner_options(cmdline::parser& cl) {
  cl.add<cmdline::path>("config",
                        'c',
                        "load config file",
                        cmdline::dontsave,
                        cmdline::path(),
                        cmdline::default_reader<cmdline::path>(),
                        cmdline::load);

  cl.add<double>(
      "r-distance", 'r', "R distance between scintillator", false, 0.5);
  cl.add<double>("s-length", 'l', "scintillator length", false, 1);
  cl.add<double>("s-pixel", 'p', "pixel size", false, 0.01);
  cl.add<int>("n-pixels", 'n', "number of pixels", cmdline::dontsave, 0);
  cl.add<int>("n-z-pixels", 0, "number of z pixels", false);
  cl.add<int>("n-y-pixels", 0, "number of y pixels", false);
  cl.add<double>(
      "s-z", 0, "TOF sigma along z axis", cmdline::alwayssave, 0.015);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
}

static void add_common_options(cmdline::parser& cl) {
  cl.add<int>("n-emissions",
              'e',
              "number of emissions",
              false,
              0,
              cmdline::not_from_file);
  cl.add("verbose", 'v', "print progress information (-v) or benchmark (-vv)");

  Common::add_openmp_options(cl);
}

void add_reconstruction_options(cmdline::parser& cl) {
  add_scanner_options(cl);
  add_common_options(cl);

  std::ostringstream msg;
  msg << "response ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in milimeters.";
  cl.footer(msg.str());

  cl.add("dl-is-time", 0, "delta is time, convert using speed of light");
  cl.add<double>("speed-of-light", 0, "speed of light", false, 299792458.);
  cl.add<double>("scale-length", 0, "scale input length values", false, 1.);
  cl.add<double>("scale-time", 0, "scale input time values", false, 1.);
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);
  cl.add<cmdline::path>("output",
                        'o',
                        "output files prefix (png)",
                        false,
                        cmdline::path(),
                        cmdline::not_from_file);
  cl.add<double>("png-max", 0, "maximum value mapped to 255 in PNG", false, 0);

  Common::add_cuda_options(cl);
}

void add_phantom_options(cmdline::parser& cl) {
  add_scanner_options(cl);
  add_common_options(cl);
  cl.add<double>(
      "scale", 0, "Scale phantom with given constant", cmdline::alwayssave, 1);

  cl.footer("phantom_description");
  cl.add<cmdline::path>("output",
                        'o',
                        "responses file",
                        false,
                        "phantom.bin",
                        cmdline::not_from_file);
}

void calculate_scanner_options(cmdline::parser& parser, int) {
  if (parser.exist("n-pixels")) {
    parser.get<int>("n-z-pixels") = parser.get<int>("n-pixels");
    parser.get<int>("n-y-pixels") = parser.get<int>("n-pixels");
  } else {
    auto R_distance = parser.get<double>("r-distance");
    auto scintillator_length = parser.get<double>("s-length");
    auto s_pixel = parser.get<double>("s-pixel");
    if (!parser.exist("n-z-pixels"))
      parser.get<int>("n-z-pixels") = std::ceil(scintillator_length / s_pixel);
    if (!parser.exist("n-y-pixels"))
      parser.get<int>("n-y-pixels") = std::ceil(2 * R_distance / s_pixel);
  }
}

}  // Strip
}  // PET2D

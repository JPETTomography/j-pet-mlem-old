#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/variant.h"

#include "options.h"

void set_detector_options(cmdline::parser& parser) {

  parser.add<double>(
      "r-distance", 'r', "R distance between scientilators", false, 500);
  parser.add<double>("s-length", 'l', "scentilator length", false, 1000);
  parser.add<double>("p-size", 'p', "pixel size", false, 5);
  parser.add<int>("n-pixels", 'n', "number of pixels", false, 200);
  parser.add<int>("n-z-pixels", '\0', "number of z pixels", false);
  parser.add<int>("n-y-pixels", '\0', "number of y pixels", false);
  parser.add<double>("s-z", 's', "Sigma z error", false, 10);
  parser.add<double>("s-dl", 'd', "Sigma dl error", false, 42);
}

void set_options_for_reconstruction(cmdline::parser& parser) {
  std::ostringstream msg;
  msg << "events_file ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in milimeters.";
  parser.footer(msg.str());

  set_detector_options(parser);

  parser.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  parser.add<int>(
      "iterations", 'I', "number of iterations (per block)", false, 1);
  parser.add<cmdline::path>(
      "output", 'o', "output files prefix (png)", false, "rec");

  parser.add("verbose", 'v', "prints the iterations information on std::out");
#if HAVE_CUDA
  parser.add("gpu", 'G', "run on GPU (via CUDA)");
  parser.add<int>("cuda-device", 'D', "CUDA device", cmdline::dontsave, 0);
  parser.add<int>("cuda-blocks", 'B', "CUDA blocks", cmdline::dontsave, 32);
  parser.add<int>(
      "cuda-threads", 'W', "CUDA threads per block", cmdline::dontsave, 512);
#endif
#if _OPENMP
  parser.add<int>(
      "n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#endif
}

void set_options_for_phantom(cmdline::parser& parser) {
  parser.footer("phantom_description");

  parser.add<cmdline::path>(
      "output", 'o', "output events file", false, "phantom.bin");
  set_detector_options(parser);
  parser.add<double>("emissions", 'e', "number of emissions", false, 500000);
#if _OPENMP
  parser.add<int>("n-threads", 'T', "number of OpenMP threads", false, 4);
}
#endif

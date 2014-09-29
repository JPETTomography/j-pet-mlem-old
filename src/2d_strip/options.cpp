#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/variant.h"

#include"options.h"

void set_options(cmdline::parser& parser) {
  std::ostringstream msg;
  msg << "events_file ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in milimeters.";
  parser.footer(msg.str());

  parser.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  parser.add<int>(
      "iterations", 'I', "number of iterations (per block)", false, 1);
  parser.add<cmdline::path>(
      "output", 'o', "output files prefix (png)", false, "rec");
  parser.add<double>(
      "r-distance", 'r', "R distance between scientilators", false, 500);
  parser.add<double>("s-length", 'l', "scentilator length", false, 1000);
  parser.add<double>("p-size", 'p', "pixel size", false, 5);
  parser.add<int>("n-pixels", 'n', "number of pixels", false, 200);
  parser.add<double>("s-z", 's', "Sigma z error", false, 10);
  parser.add<double>("s-dl", 'd', "Sigma dl error", false, 42);
  parser.add<double>("gm", 'u', "Gamma error", false, 0);
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

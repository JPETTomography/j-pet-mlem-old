#include <iostream>
#include <vector>
#include <ctime>
#include <sstream>
#include <iomanip>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/png_writer.h"
#include "util/progress.h"
#include "util/variant.h"

#include "event.h"
#include "reconstruction.h"

#if HAVE_CUDA
#include "cuda/reconstruction.h"
#endif

using namespace std;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    std::ostringstream msg;
    msg << "events_file ..." << std::endl;
    msg << "build: " << VARIANT << std::endl;
    msg << "note: All length options below should be expressed in meters.";
    cl.footer(msg.str());

    cl.add<int>("i-blocks", 'i', "number of iteration blocks", false, 0);
    cl.add<int>(
        "iterations", 'I', "number of iterations (per block)", false, 1);
    cl.add<cmdline::path>(
        "output", 'o', "output files prefix (png)", false, "rec");
    cl.add<double>(
        "r-distance", 'r', "R distance between scientilators", false, 500);
    cl.add<double>("s-length", 'l', "scentilator length", false, 1000);
    cl.add<double>("p-size", 'p', "pixel size", false, 5);
    cl.add<int>("n-pixels", 'n', "number of pixels", false, 200);
    cl.add<double>("s-z", 's', "Sigma z error", false, 10);
    cl.add<double>("s-dl", 'd', "Sigma dl error", false, 63);
    cl.add<double>("gm", 'u', "Gamma error", false, 0);
#if HAVE_CUDA
    cl.add("gpu", 'G', "run on GPU (via CUDA)");
    cl.add<int>("cuda-device", 'D', "CUDA device", cmdline::dontsave, 0);
    cl.add<int>("cuda-blocks", 'B', "CUDA blocks", cmdline::dontsave, 1);
    cl.add<int>(
        "cuda-threads", 'W', "CUDA threads per block", cmdline::dontsave, 512);
#endif
#if _OPENMP
    cl.add<int>(
        "n-threads", 't', "number of OpenMP threads", cmdline::dontsave);
#endif

    cl.parse_check(argc, argv);

    if (!cl.rest().size()) {
      throw("at least one events input file expected, consult --help");
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    auto R_distance = cl.get<double>("r-distance");
    auto scintillator_length = cl.get<double>("s-length");
    auto pixel_size = cl.get<double>("p-size");
    auto n_pixels = scintillator_length / pixel_size;
    auto sigma = cl.get<double>("s-z");
    auto dl = cl.get<double>("s-dl");

    Reconstruction<double> reconstruction(
        R_distance, scintillator_length, n_pixels, pixel_size, sigma, dl);

    for (auto& fn : cl.rest()) {
      ibstream in(fn);
      reconstruction << in;
    }

    auto n_blocks = cl.get<int>("i-blocks");
    auto n_iterations = cl.get<int>("iterations");
    auto output_wo_ext = cl.get<cmdline::path>("output").wo_ext();

    Progress progress(true, n_blocks * n_iterations, 1);

#if HAVE_CUDA
    if (cl.exist("gpu")) {
      Reconstruction<float> sp_reconstruction(
          R_distance, scintillator_length, n_pixels, pixel_size, sigma, dl);
      run_gpu_reconstruction(sp_reconstruction.detector,
                             reconstruction.get_event_list(),
                             n_blocks,
                             n_iterations,
                             cl.get<int>("cuda-device"),
                             cl.get<int>("cuda-blocks"),
                             cl.get<int>("cuda-threads"),
                             progress,
                             output_wo_ext);
    } else
#endif
    {
      for (int block = 0; block < n_blocks; block++) {
        reconstruction(progress, n_iterations, block * n_iterations);

        std::stringstream fn;
        fn << output_wo_ext << "_"               // phantom_
           << std::setw(3) << std::setfill('0')  //
           << block * n_iterations + 1           // 001
           << std::setw(0) << ".png";            // .png

        png_writer png(fn.str());
        reconstruction.output_bitmap(png);
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}

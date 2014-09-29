#include <iostream>
#include <ostream>

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

#include "options.h"
#include "event.h"
#include "reconstruction.h"

#if HAVE_CUDA
#include "cuda/reconstruction.h"
#endif


std::ostream& print_statistics(std::ostream& out,
                               const Reconstruction<double>& reconstruction,
                               int n_iterations,
                               int n_blocks);

using namespace std;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    set_options_for_reconstruction(cl);
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

    auto sigma = cl.get<double>("s-z");
    auto dl = cl.get<double>("s-dl");

    int n_z_pixels = (int)std::ceil(scintillator_length / pixel_size);
    int n_y_pixels = (int)std::ceil(2 * R_distance / pixel_size);
    std::cerr << n_y_pixels << "x" << n_z_pixels << std::endl;
    Reconstruction<double> reconstruction(R_distance,
                                          scintillator_length,
                                          n_y_pixels,
                                          n_z_pixels,
                                          pixel_size,
                                          pixel_size,
                                          sigma,
                                          dl);

    for (auto& fn : cl.rest()) {
      ibstream in(fn);
      reconstruction << in;
    }

    auto n_blocks = cl.get<int>("blocks");
    auto n_iterations = cl.get<int>("iterations");
    auto output_wo_ext = cl.get<cmdline::path>("output").wo_ext();

    Progress progress(true, n_blocks * n_iterations, 1);

#if HAVE_CUDA
    if (cl.exist("gpu")) {
      StripDetector<float> detector(R_distance,
                                    scintillator_length,
                                    n_y_pixels,
                                    n_z_pixels,
                                    pixel_size,
                                    pixel_size,
                                    sigma,
                                    dl);
      run_gpu_reconstruction(detector,
                             reconstruction.get_event_list(),
                             n_blocks,
                             n_iterations,
                             cl.get<int>("cuda-device"),
                             cl.get<int>("cuda-blocks"),
                             cl.get<int>("cuda-threads"),
                             cl.exist("verbose"),
                             progress,
                             output_wo_ext);
    } else
#endif
    {
      png_writer png(output_wo_ext + "_sensitivity.png");
      reconstruction.output_bitmap(png, true);

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

      if (cl.exist("verbose")) {
        print_statistics(std::cout, reconstruction, n_iterations, n_blocks);
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}

std::ostream& print_statistics(std::ostream& out,
                               const Reconstruction<double>& reconstruction,
                               int n_iterations,
                               int n_blocks) {
  size_t iterations = n_iterations * n_blocks;
  size_t events = reconstruction.n_events_processed() / iterations;
  size_t pixels = reconstruction.n_pixels_processed() / iterations;
  size_t kernels = reconstruction.n_kernel_calls() / iterations;
  out << "iterations: " << iterations << " "
      << "events: " << events << " "
      << "pixels: " << pixels << " "
      << "(" << (double)pixels / events << ") "
      << "kernel calls: " << kernels << " "
      << "(" << (double)kernels / events << ")" << std::endl;

  size_t bb_width_sum = reconstruction.bb_width_sum() / n_iterations;
  size_t bb_height_sum = reconstruction.bb_height_sum() / n_iterations;
  size_t bb_width2_sum = reconstruction.bb_width2_sum() / n_iterations;
  size_t bb_height2_sum = reconstruction.bb_height2_sum() / n_iterations;
  size_t bb_width_height_sum =
      reconstruction.bb_width_height_sum() / n_iterations;
  double avg_width = (double)bb_width_sum / events;
  double avg_height = (double)bb_height_sum / events;
  double avg_width2 = (double)bb_width2_sum / events;
  double avg_height2 = (double)bb_height2_sum / events;
  double avg_width_height = (double)bb_width_height_sum / events;
  avg_width2 -= avg_width * avg_width;
  avg_height2 -= avg_height * avg_height;
  avg_width_height -= avg_width * avg_height;
  out << "width: " << avg_width << "(" << std::sqrt(avg_width2) << ")"
      << " height: " << avg_height << "(" << std::sqrt(avg_height2) << ")  "
      << avg_width_height << std::endl;
}

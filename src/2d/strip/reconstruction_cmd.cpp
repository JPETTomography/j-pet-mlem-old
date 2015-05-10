/// \page cmd_2d_strip_reconstruction 2d_strip_reconstruction
/// \brief 2D Strip PET reconstruction tool
///
/// Reconstructs image using List-Mode with analytic kernel approximation from
/// physical scanner response or simulated response output from \ref
/// cmd_2d_strip_phantom.
///
/// \image html cs000_ev.pdf.png
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/strip/reconstruction_cmd.txt
///
/// \sa \ref cmd_2d_strip_phantom

#include <iostream>
#include <ostream>

#include <vector>
#include <sstream>
#include <iomanip>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/png_writer.h"
#include "util/progress.h"

#include "options.h"
#include "event.h"
#include "reconstruction.h"

#if HAVE_CUDA
#include "cuda/reconstruction.h"
#endif

using namespace PET2D;
using namespace PET2D::Strip;

void print_statistics(std::ostream& out,
                      const Reconstruction<double>& reconstruction,
                      int n_iterations,
                      int n_blocks,
                      std::string prefix = std::string());

using namespace std;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    add_reconstruction_options(cl);
    cl.parse_check(argc, argv);
    calculate_scanner_options(cl);

    if (!cl.rest().size()) {
      throw("at least one events input file expected, consult --help");
    }

    cmdline::load_accompanying_config(cl);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    Scanner<float, short> scanner(PET2D_STRIP_SCANNER_CL(cl));
    Reconstruction<float> reconstruction(scanner);

    auto verbose = cl.exist("verbose");
    auto verbosity = cl.count("verbose");
    if (verbose) {
      std::cout << "# image: " << scanner.n_y_pixels << "x"
                << scanner.n_z_pixels << std::endl;
    }

    for (auto& fn : cl.rest()) {
      if (cmdline::path(fn).ext() == ".txt") {
        std::ifstream events(fn);
        if (!events.is_open()) {
          throw("cannot open phantom events file: " + fn);
        }
        reconstruction << events;
      } else {
        util::ibstream events(fn);
        if (!events.is_open()) {
          throw("cannot open phantom events file: " + fn);
        }
        reconstruction << events;
      }
      if (verbose) {
        std::cerr << "# read " << reconstruction.events.size()
                  << " events from " << fn << std::endl;
      }

      auto dl_is_time = cl.exist("dl-is-time");
      if (cl.exist("scale-length")) {
        auto scale = cl.get<double>("scale-length");
        for (auto& event : reconstruction.events) {
          event.z_u *= scale;
          event.z_d *= scale;
          if (!dl_is_time)
            event.dl *= scale;
        }
      }
      if (dl_is_time) {
        auto scale = cl.get<double>("scale-time");
        auto speed_of_light = cl.get<double>("speed-of-light");
        for (auto& event : reconstruction.events) {
          event.dl = event.dl * speed_of_light * scale;
        }
      }
    }

    auto n_blocks = cl.get<int>("blocks");
    auto n_iterations = cl.get<int>("iterations");
    auto output_name = cl.get<cmdline::path>("output");
    auto output_base_name = output_name.wo_ext();
    auto output_txt = output_name.ext() == ".txt";

    util::progress progress(verbosity, n_blocks * n_iterations, 1);

    int n_file_digits = n_blocks * n_iterations >= 1000 ? 4 :
                        n_blocks * n_iterations >= 100 ? 3 :
                        n_blocks * n_iterations >= 10 ? 2 : 1;

#if HAVE_CUDA
    if (cl.exist("gpu")) {
      Scanner<float, short> single_precision_scanner(scanner);

      GPU::run_reconstruction(single_precision_scanner,
                              reconstruction.events,
                              n_blocks,
                              n_iterations,
                              cl.get<int>("cuda-device"),
                              cl.get<int>("cuda-blocks"),
                              cl.get<int>("cuda-threads"),
                              cl.exist("verbose"),
                              progress,
                              output_base_name,
                              output_txt,
                              n_file_digits);
    } else
#endif
    {
      if (output_base_name.length()) {
        util::png_writer png(output_base_name + "_sensitivity.png");
        reconstruction.output_bitmap(png, true);
      }

      for (int block = 0; block < n_blocks; block++) {
        reconstruction(progress, n_iterations, block * n_iterations);

        if (output_base_name.length()) {
          std::stringstream fn;
          fn << output_base_name << "_"      // phantom_
             << std::setw(n_file_digits)     //
             << std::setfill('0')            //
             << (block + 1) * n_iterations;  // 001

          util::png_writer png(fn.str() + ".png");
          reconstruction.output_bitmap(png);

          if (output_txt) {
            std::ofstream txt(fn.str() + ".txt");
            txt.precision(12);
            txt << std::fixed;
            reconstruction.output_tuples(txt);
          } else {
            util::obstream bin(fn.str() + ".bin");
            reconstruction >> bin;
          }
        }
      }
#if USE_STATISTICS
      if (verbose) {
        print_statistics(
            std::cout, reconstruction, n_iterations, n_blocks, "# ");
      }
#endif
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}

#if USE_STATISTICS
void print_statistics(std::ostream& out,
                      const Reconstruction<double>& reconstruction,
                      int n_iterations,
                      int n_blocks,
                      std::string prefix) {
  size_t iterations = n_iterations * n_blocks;
  size_t events = reconstruction.stats.n_events_processed / iterations;
  size_t pixels = reconstruction.stats.n_pixels_processed / iterations;
  size_t kernels = reconstruction.stats.n_kernel_calls / iterations;
  out << prefix << "iterations: " << iterations << " "
      << "events: " << events << " "
      << "pixels: " << pixels << " "
      << "(" << (double)pixels / events << ") "
      << "kernel calls: " << kernels << " "
      << "(" << (double)kernels / events << ")" << std::endl;

  size_t bb_width_sum = reconstruction.stats.bb_width_sum / n_iterations;
  size_t bb_height_sum = reconstruction.stats.bb_height_sum / n_iterations;
  size_t bb_width2_sum = reconstruction.stats.bb_width2_sum / n_iterations;
  size_t bb_height2_sum = reconstruction.stats.bb_height2_sum / n_iterations;
  size_t bb_width_height_sum =
      reconstruction.stats.bb_width_height_sum / n_iterations;
  double avg_width = (double)bb_width_sum / events;
  double avg_height = (double)bb_height_sum / events;
  double avg_width2 = (double)bb_width2_sum / events;
  double avg_height2 = (double)bb_height2_sum / events;
  double avg_width_height = (double)bb_width_height_sum / events;
  avg_width2 -= avg_width * avg_width;
  avg_height2 -= avg_height * avg_height;
  avg_width_height -= avg_width * avg_height;
  out << prefix << "width: " << avg_width << "(" << std::sqrt(avg_width2) << ")"
      << " height: " << avg_height << "(" << std::sqrt(avg_height2) << ")  "
      << avg_width_height << std::endl;
}
#endif

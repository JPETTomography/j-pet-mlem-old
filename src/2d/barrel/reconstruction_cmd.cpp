/// \page cmd_2d_barrel_reconstruction 2d_barrel_reconstruction
/// \brief 2D Barrel PET reconstruction tool
///
/// Reconstructs image using given system matrix produced by \ref
/// cmd_2d_barrel_matrix and mean file representing physical detector response
/// or
/// simulated response output from \ref cmd_2d_barrel_phantom.
///
/// Authors
/// -------
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
///
/// References
/// ----------
/// Based on:
///  "Implementing and Accelerating the EM Algorithm for Positron Emission
///   Tomography" by Linda Kaufman
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/reconstruction_cmd.txt
///
/// \sa \ref cmd_2d_barrel_matrix, \ref cmd_2d_barrel_phantom

#ifdef __SSE3__
#include <xmmintrin.h>
#include <pmmintrin.h>
#endif

#include <ctime>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "reconstruction.h"
#include "util/png_writer.h"
#include "util/backtrace.h"
#include "options.h"

#if _OPENMP
#include <omp.h>
#endif

using F = float;
using S = short;
using Hit = int;

using Reconstruction = PET2D::Barrel::Reconstruction<F, S, Hit>;

int main(int argc, char* argv[]) {

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    PET2D::Barrel::add_reconstruction_options(cl);
    cl.parse_check(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    auto verbose = cl.exist("verbose");
    auto use_sensitivity = !cl.exist("no-sensitivity");

    util::ibstream in_matrix(cl.get<cmdline::path>("system"));
    if (!in_matrix.is_open())
      throw("cannot open input file: " + cl.get<cmdline::path>("system"));
    Reconstruction::Matrix matrix(in_matrix);

    if (matrix.triangular()) {
      throw("matrix must be in full form");
    }

    std::ifstream in_means(cl.get<cmdline::path>("mean"));
    if (!in_means.is_open())
      throw("cannot open input file: " + cl.get<cmdline::path>("mean"));

    int n_i_blocks = cl.get<int>("blocks");
    Reconstruction reconstruction(matrix, in_means, use_sensitivity);

    if (verbose) {
      std::cerr << "reconstruction:" << std::endl;
#if _OPENMP
      std::cerr << " threads       = " << omp_get_max_threads() << std::endl;
#endif
      std::cerr << " pixels in row = " << matrix.n_pixels_in_row() << std::endl;
      std::cerr << " TOF positions = " << matrix.n_tof_positions() << std::endl;
      std::cerr << " emissions     = " << matrix.n_emissions() << std::endl;
      std::cerr << " detectors     = " << matrix.n_detectors() << std::endl;
    }

    auto n_pixels_in_row = reconstruction.n_pixels_in_row();
    auto total_n_pixels = n_pixels_in_row * n_pixels_in_row;
    Reconstruction::Output rho(total_n_pixels, 0.0);
    Reconstruction::Output rho_detected(total_n_pixels, 0.0);

    // no output, just make reconstruction in place and exit!
    if (!cl.exist("output")) {
      for (int i = 0; i < n_i_blocks; ++i) {
        reconstruction.emt(cl.get<int>("iterations"));
        rho = reconstruction.rho();
      }
      return 0;
    }

    auto output = cl.get<cmdline::path>("output");

    std::ofstream out;
    std::ofstream out_detected;

    out.open(output);
    out_detected.open(output.wo_ext() + "_detected" + output.ext());

    double sec = 0.0;
    auto n_iterations = cl.get<int>("iterations");

    for (int i = 0; i < n_i_blocks; ++i) {
      clock_t start = clock();
      reconstruction.emt(n_iterations);
      clock_t stop = clock();
      sec += static_cast<double>(stop - start) / CLOCKS_PER_SEC;

      out << reconstruction.rho();
      out_detected << reconstruction.rho_detected();
    }

    if (use_sensitivity) {
      std::ofstream out_sensitivity(output.wo_ext() + "_sensitivity" +
                                    output.ext());
      int n_row = 0;
      for (auto& scale : reconstruction.scale()) {
        if (scale > 0)
          out_sensitivity << 1 / scale;
        else
          out_sensitivity << 0;

        if (++n_row >= n_pixels_in_row) {
          n_row = 0;
          out_sensitivity << "\n";
        } else {
          out_sensitivity << " ";
        }
      }

      util::png_writer png(output.wo_ext() + "_sensitivity.png");
      png << reconstruction.sensitivity();
    }

    std::cout << std::endl;
    if (verbose) {
      std::cout << "time = " << sec << "s "
                << "time/iter = " << sec / (n_i_blocks * n_iterations) << "s"
                << std::endl;
    }

    out.close();
    out_detected.close();

    // output reconstruction PNG
    util::png_writer png(output.wo_ext() + ".png");
    png << rho;

    return 0;
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  }
  return 1;
}

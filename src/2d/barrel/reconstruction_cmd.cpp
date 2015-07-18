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

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "reconstruction.h"
#include "util/png_writer.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "options.h"

#if _OPENMP
#include <omp.h>
#endif

using F = float;
using S = short;
using Hit = int;

using Reconstruction = PET2D::Barrel::Reconstruction<F, S, Hit>;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  cmdline::parser cl;
  PET2D::Barrel::add_reconstruction_options(cl);
  cl.parse_check(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto verbose = cl.count("verbose");
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

  auto n_blocks = cl.get<int>("blocks");
  auto n_iterations = cl.get<int>("iterations");

  Reconstruction reconstruction(matrix, in_means, use_sensitivity);
  auto n_pixels_in_row = reconstruction.n_pixels_in_row();

  // no output, just make reconstruction in place and exit!
  if (!cl.exist("output")) {
    util::progress progress(verbose, n_blocks * n_iterations, 1);
    for (int block = 0; block < n_blocks; ++block) {
      reconstruction(progress, n_iterations, block * n_iterations);
    }
    return 0;
  }

  auto output = cl.get<cmdline::path>("output");
  std::ofstream out_rho(output);

  util::progress progress(verbose, n_blocks * n_iterations, 1);
  for (int block = 0; block < n_blocks; ++block) {
    reconstruction(progress, n_iterations, block * n_iterations);
    out_rho << reconstruction.rho();
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

  // output reconstruction PNG
  util::png_writer png(output.wo_ext() + ".png");
  png << reconstruction.rho();

  CMDLINE_CATCH
}

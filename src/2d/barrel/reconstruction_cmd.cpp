/// \page cmd_2d_barrel_reconstruction 2d_barrel_reconstruction
/// \brief 2D Barrel PET reconstruction tool
///
/// Reconstructs image using given system matrix produced by \ref
/// cmd_2d_barrel_matrix and mean file representing physical detector response
/// or simulated response output from \ref cmd_2d_barrel_phantom.
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
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"
#include "util/cmdline_types.h"
#include "reconstruction.h"
#include "util/png_writer.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "options.h"

#include "common/types.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/reconstruction.h"
#include "simple_geometry.h"
#endif

using Reconstruction = PET2D::Barrel::Reconstruction<F, S, Hit>;
#if HAVE_CUDA
using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
#endif

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#ifdef __SSE3__
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif

  cmdline::parser cl;
  cl.footer("means");
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

  Reconstruction reconstruction(matrix, use_sensitivity);

  for (const auto& fn : cl.rest()) {
    std::ifstream in_means(fn);
    if (!in_means.is_open())
      throw("cannot open input file: " + cl.get<cmdline::path>("mean"));
    reconstruction << in_means;
  }
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
  auto output_base_name = output.wo_ext();

  util::progress progress(verbose, n_blocks * n_iterations, 1);

  const int n_file_digits = n_blocks * n_iterations >= 1000
                                ? 4
                                : n_blocks * n_iterations >= 100
                                      ? 3
                                      : n_blocks * n_iterations >= 10 ? 2 : 1;

#if HAVE_CUDA
  if (cl.exist("gpu")) {
    SimpleGeometry geometry(matrix);
    PET2D::Barrel::GPU::Reconstruction::run(
        geometry,
        reconstruction.means().data(),
        reconstruction.means().size(),
        n_pixels_in_row,
        n_pixels_in_row,
        n_blocks,
        n_iterations,
        [&](int iteration, float* output) {
          std::stringstream fn;
          fn << output_base_name << "_";  // phantom_
          if (iteration >= 0) {
            fn << std::setw(n_file_digits)    //
               << std::setfill('0')           //
               << iteration << std::setw(0);  // 001
          } else {
            fn << "sensitivity";
          }

          util::png_writer png(fn.str() + ".png");
          png.write(n_pixels_in_row, n_pixels_in_row, output);
        },
        [&](int completed, bool finished) { progress(completed, finished); },
        cl.get<int>("cuda-device"),
        cl.get<int>("cuda-blocks"),
        cl.get<int>("cuda-threads"),
        [](const char* device_name) {
          std::cerr << "   CUDA device = " << device_name << std::endl;
        });
  } else
#endif
  {
    std::ofstream out_rho(output);
    for (int block = 0; block < n_blocks; ++block) {
      reconstruction(progress, n_iterations, block * n_iterations);
      out_rho << reconstruction.rho();
    }
  }

  if (use_sensitivity) {
    std::ofstream out_sensitivity(output_base_name + "_sensitivity" +
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

    util::png_writer png(output_base_name + "_sensitivity.png");
    png << reconstruction.sensitivity();
  }

  // output reconstruction PNG
  util::png_writer png(output_base_name + ".png", cl.get<double>("png-max"));
  png << reconstruction.rho();

  util::png_writer png_detected(output_base_name + "_detected.png",
                                cl.get<double>("png-max"));
  png_detected << reconstruction.rho_detected();

  CMDLINE_CATCH
}

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
/// \verboutput 2d_barrel_reconstruction
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
#include "geometry_soa.h"
#endif

using Reconstruction = PET2D::Barrel::Reconstruction<F, S, Hit>;
#if HAVE_CUDA
using GeometrySOA = PET2D::Barrel::GeometrySOA<F, S>;
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
  ENSURE_IS_OPEN(in_matrix, "input matrix", cl.get<cmdline::path>("system"));
  Reconstruction::Matrix matrix(in_matrix);

  if (matrix.triangular()) {
    throw(
        "matrix must be in full form, "
        "convert using 2d_barrel_matrix --full option");
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
  auto n_iterations_in_block = cl.get<int>("iterations");
  auto n_iterations = n_blocks * n_iterations_in_block;

  Reconstruction reconstruction(matrix, use_sensitivity);
  auto n_pixels_in_row = reconstruction.n_pixels_in_row();

  for (const auto& fn : cl.rest()) {
    std::ifstream in_means(fn);
    ENSURE_IS_OPEN(in_means, "means", fn);
    reconstruction << in_means;
  }

  auto output_name = cl.get<cmdline::path>("output");
  auto output_base_name = output_name.wo_ext();
  auto output_ext = output_name.ext();
  auto output_txt = output_ext == ".txt";

  util::progress progress(verbose, n_iterations, 1);

#if HAVE_CUDA
  if (cl.exist("gpu")) {
    GeometrySOA geometry(matrix);
    PET2D::Barrel::GPU::Reconstruction::run(
        geometry,
        reconstruction.means().data(),
        reconstruction.means().size(),
        n_pixels_in_row,
        n_pixels_in_row,
        n_blocks,
        n_iterations_in_block,
        [&](int iteration,
            const PET2D::Barrel::GPU::Reconstruction::Output& output) {
          if (!output_base_name.length())
            return;
          auto fn = iteration >= 0
                        ? output_base_name.add_index(iteration, n_iterations)
                        : output_base_name + "_sensitivity";
          util::png_writer png(fn + ".png");
          png << output;
          if (output_txt) {
            std::ofstream txt(fn + ".txt");
            txt << output;
          } else if (output_ext != ".png") {
            util::obstream bin(fn + output_ext);
            bin << output;
          }
        },
        [&](int completed, bool finished) { progress(completed, finished); },
        cl.get<int>("cuda-device"),
        cl.get<int>("cuda-blocks"),
        cl.get<int>("cuda-threads"),
        [=](const char* device_name) {
          if (verbose) {
            std::cerr << "   CUDA device = " << device_name << std::endl;
          }
        });
  } else
#endif
  {
    if (use_sensitivity && output_base_name.length()) {
      std::ofstream out_sensitivity(output_base_name + "_sensitivity" +
                                    output_ext);
      // FIXME: simplify sensitivity output here
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

    for (int block = 0; block < n_blocks; ++block) {
      reconstruction(
          progress, n_iterations_in_block, block * n_iterations_in_block);
      if (!output_base_name.length())
        continue;
      auto fn = output_base_name.add_index((block + 1) * n_iterations_in_block,
                                           n_iterations);
      util::png_writer png(fn + ".png");
      png << reconstruction.rho();
      if (output_txt) {
        std::ofstream txt(fn + ".txt");
        txt << reconstruction.rho();
      } else if (output_ext != ".png") {
        util::obstream bin(fn + output_ext);
        bin << reconstruction.rho();
      }
    }
  }

  CMDLINE_CATCH
}

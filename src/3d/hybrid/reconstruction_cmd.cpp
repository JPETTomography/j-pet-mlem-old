/// \page cmd_3d_hybrid_reconstruction 3d_hybrid_reconstruction
/// \brief 3D Longitudinal PET reconstruction tool
///
/// Reconstructs image using given system matrix produced by \ref
/// cmd_3d_hybrid_matrix and mean file representing physical detector response
/// or simulated response output from \ref cmd_3d_hybrid_phantom.
///
/// Usage
/// -----
/// \verbinclude src/3d/hybrid/reconstruction_cmd.txt
///
/// \sa \ref cmd_3d_hybrid_matrix, \ref cmd_3d_hybrid_phantom

#include <iostream>
#include <fstream>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/backtrace.h"

#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/lor_info.h"
#include "2d/strip/gausian_kernel.h"

#include "scanner.h"
#include "reconstruction.h"
#include "options.h"

#include "common/types.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

using RNG = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Point = PET2D::Point<F>;
using Geometry = PET2D::Barrel::Geometry<F, S>;
using MathematicaGraphics = Common::MathematicaGraphics<F>;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Hybrid::add_reconstruction_options(cl);
  cl.parse_check(argc, argv);
  PET3D::Hybrid::calculate_scanner_options(cl, argc);

  Scanner scanner(
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, F)),
      F(cl.get<double>("length")));
  scanner.set_sigmas(cl.get<double>("s-z"), cl.get<double>("s-dl"));
  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto lor_info_file_name = cl.get<std::string>("geometry");
  util::ibstream in_geometry(lor_info_file_name);
  Geometry geometry(in_geometry);

  if (geometry.n_detectors != (int)scanner.barrel.size()) {
    throw("n_detectors mismatch");
  }

  if (cl.exist("verbose")) {
    std::cout << geometry.n_detectors << std::endl;
    std::cout << geometry.grid.n_columns << "x" << geometry.grid.n_rows << " "
              << geometry.grid.pixel_size << std::endl;
  }

  PET3D::Hybrid::Reconstruction<Scanner, PET2D::Strip::GaussianKernel<F>>
      reconstruction(scanner, geometry, -0.200, 80);

  reconstruction.calculate_weight();
  reconstruction.calculate_sensitivity();

  for (const auto& fn : cl.rest()) {
    std::ifstream in_response(fn);
    reconstruction << in_response;
  }

  {
    std::ofstream out_graphics("event.m");
    MathematicaGraphics graphics(out_graphics);
    graphics.add(scanner.barrel);
    reconstruction.graph_frame_event(graphics, 0);
  }

  auto n_blocks = cl.get<int>("blocks");
  auto n_iter = cl.get<int>("iterations");

  for (int block = 0; block < n_blocks; ++block) {
    for (int i = 0; i < n_iter; i++) {
      std::cout << block * n_iter + i << " " << reconstruction() << "\n";
    }
    char rho_file_name[64];
    sprintf(rho_file_name,
            "%s_%03d.bin",
            output_base_name.c_str(),
            (block + 1) * n_iter);
    util::obstream out(rho_file_name);
    out << reconstruction.rho();
  }
  std::cout << reconstruction.event_count() << " "
            << reconstruction.voxel_count() << " "
            << reconstruction.pixel_count() << "\n";
  std::cout << (double)reconstruction.voxel_count() /
                   reconstruction.event_count()
            << " ";
  std::cout << (double)reconstruction.pixel_count() /
                   reconstruction.event_count()
            << "\n";

  util::obstream out(output);
  out << reconstruction.rho();

  CMDLINE_CATCH
}

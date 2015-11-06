/// \page cmd_3d_hybrid_reconstruction 3d_hybrid_reconstruction
/// \brief 3D Longitudinal PET reconstruction tool
///
/// Reconstructs image using given system matrix produced by \ref
/// cmd_3d_hybrid_matrix and mean file representing physical detector response
/// or simulated response output from \ref cmd_3d_hybrid_phantom.
///
/// Usage
/// -----
/// \verboutput 3d_hybrid_reconstruction
///
/// \sa \ref cmd_3d_hybrid_matrix, \ref cmd_3d_hybrid_phantom

#include <iostream>
#include <fstream>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/backtrace.h"
#include "util/png_writer.h"
#include "util/nrrd_writer.h"
#include "util/progress.h"

#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/lor_geometry.h"
#include "2d/strip/gausian_kernel.h"

#include "2d/barrel/geometry_matrix_loader.h"
#include "2d/barrel/sparse_matrix.h"

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

#if HAVE_CUDA
#include "cuda/reconstruction.h"
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
  cl.add<float>("z-left", 0, "left extent in z direction");
  cl.add<int>("n-planes", 0, "number of voxels in z direction");
  cl.add<cmdline::path>("system", 0, "system matrix file", false);
  cl.parse_check(argc, argv);
  PET3D::Hybrid::calculate_scanner_options(cl, argc);

  auto verbose = cl.count("verbose");

  Scanner scanner(
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, F)),
      F(cl.get<double>("length")));
  scanner.set_sigmas(cl.get<double>("s-z"), cl.get<double>("s-dl"));

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  util::ibstream in_geometry(cl.get<std::string>("geometry"));
  Geometry geometry(in_geometry);

  if (geometry.n_detectors != (int)scanner.barrel.size()) {
    throw("n_detectors mismatch");
  }

  if (cl.exist("verbose")) {
    std::cout << geometry.n_detectors << std::endl;
    std::cout << geometry.grid.n_columns << "x" << geometry.grid.n_rows << " "
              << geometry.grid.pixel_size << std::endl;
  }

  if (cl.exist("system")) {
    load_system_matrix_from_file<F, S, Hit>(
        cl.get<cmdline::path>("system"), geometry, verbose);
    std::cerr << "loaded system matrix" << std::endl;
  }

  PET3D::Hybrid::Reconstruction<Scanner, PET2D::Strip::GaussianKernel<F>>
  reconstruction(
      scanner, geometry, cl.get<float>("z-left"), cl.get<int>("n-planes"));

  if (!cl.exist("system")) {
    reconstruction.calculate_weight();
  }
  reconstruction.calculate_sensitivity();
  reconstruction.normalize_geometry_weights();

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
  auto n_iterations_in_block = cl.get<int>("iterations");
  auto n_iterations = n_blocks * n_iterations_in_block;

  auto output_name = cl.get<cmdline::path>("output");
  auto output_base_name = output_name.wo_ext();
  auto output_ext = output_name.ext();
  auto output_txt = output_ext == ".txt";

  util::progress progress(verbose, n_iterations, 1);

#if HAVE_CUDA
  if (cl.exist("gpu")) {
    PET3D::Hybrid::GPU::Reconstruction::SimpleGeometry simple_geometry(
        geometry);
    PET3D::Hybrid::GPU::Reconstruction::run(
        simple_geometry,
        reconstruction.events().data(),
        reconstruction.n_events(),
        reconstruction.scanner.sigma_z(),
        reconstruction.scanner.sigma_dl(),
        geometry.grid.n_columns,
        geometry.grid.n_rows,
        reconstruction.n_planes,
        n_blocks,
        n_iterations_in_block,
        [&](int iteration,
            const PET3D::Hybrid::GPU::Reconstruction::Output& output) {
          if (!output_base_name.length())
            return;
          auto fn = iteration >= 0
                        ? output_base_name.add_index(iteration, n_iterations)
                        : output_base_name + "_sensitivity";
          util::nrrd_writer nrrd(fn + ".nrrd", fn + output_ext, output_txt);
          nrrd << output;
          if (output_txt) {
            std::ofstream txt(fn + ".txt");
            txt << output;
          } else {
            util::obstream bin(fn + output_ext);
            bin << output;
          }
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
    for (int block = 0; block < n_blocks; ++block) {
      for (int i = 0; i < n_iterations_in_block; i++) {
        std::cout << block * n_iterations_in_block + i << ' '
                  << reconstruction() << "\n";
      }
      auto fn = output_base_name.add_index((block + 1) * n_iterations_in_block,
                                           n_iterations);
      util::nrrd_writer nrrd(fn + ".nrrd", fn + output_ext, output_txt);
      nrrd << reconstruction.rho();
      if (output_txt) {
        std::ofstream txt(fn + ".txt");
        txt << reconstruction.rho();
      } else {
        util::obstream bin(fn + output_ext);
        bin << reconstruction.rho();
      }
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
  }

  CMDLINE_CATCH
}

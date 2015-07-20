/// \page cmd_2d_barrel_lm_reconstruction 2d_barrel_lm_reconstruction
/// \brief 2D Barrel PET LM reconstruction tool
///
/// Reconstructs image using LM method using given geometry description
/// produced by \ref cmd_2d_barrel_geometry and mean file representing physical
/// detector response or simulated response output from
/// \ref cmd_2d_barrel_phantom.
///
/// Authors
/// -------
/// - Piotr Bialas <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/lm_reconstruction_cmd.txt
///
/// \sa \ref cmd_2d_barrel_geometry, \ref cmd_2d_barrel_phantom

#include <iostream>
#include <fstream>
#include <random>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"
#include "util/progress.h"
#include "util/backtrace.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/options.h"
#include "2d/barrel/lor_info.h"
#include "2d/barrel/sparse_matrix.h"
#include "2d/barrel/lm_reconstruction.h"

#include "common/mathematica_graphics.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

using F = float;
using S = int16_t;
using Hit = int;
using RNG = std::mt19937;
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<Scanner>;
using MathematicaGraphics = Common::MathematicaGraphics<F>;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  cl.add<cmdline::path>("geometry", 0, "geometry information", true);
  cl.add<cmdline::path>("system", 0, "system maxtrix", false);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);

  cl.add<double>("length", 0, "length of the detector", false, 0.3);
  cl.add<cmdline::path>("response", 0, "detector responses", true);

  PET2D::Barrel::add_matrix_options(cl);
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);
  cl.add("graphics", 'g', "output mathematica .m graphics file", false);
  cl.add("event", 0, "event number", false, 0);

  cl.try_parse(argc, argv);

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto verbose = cl.count("verbose");

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  util::ibstream in_geometry(cl.get<cmdline::path>("geometry"));
  PET2D::Barrel::Geometry<F, S> geometry(in_geometry);

  if (verbose) {
    std::cout << geometry.n_detectors << std::endl;
    std::cout << geometry.grid.n_columns << "x" << geometry.grid.n_rows << " "
              << geometry.grid.pixel_size << std::endl;
  }

  if (cl.exist("system")) {
    geometry.erase_pixel_info();
    auto system_matrix_file_name = cl.get<cmdline::path>("system");
    util::ibstream system_matrix_istream(system_matrix_file_name);
    PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit> matrix(system_matrix_istream);
    if (verbose) {
      std::cout << "read in system matrix" << std::endl;
    }
    matrix.sort_by_lor_n_pixel();
    matrix.merge_duplicates();
    F n_emissions = F(matrix.n_emissions());
    if (geometry.grid.n_columns != matrix.n_pixels_in_row()) {
      throw("mismatch in number of pixels with matrix");
    }
    if (matrix.triangular()) {
      throw("matrix is not full");
    }

    for (auto& element : matrix) {
      auto lor = element.lor;
      F weight = element.hits / n_emissions;
      geometry.push_back_pixel(lor, element.pixel, weight);
    }
    geometry.sort_all();
  }

  PET2D::Barrel::LMReconstruction<F, S> reconstruction(
      geometry, cl.get<double>("s-dl") / 2);
  if (verbose)
    std::cout << "created reconstruction\n";
  if (cl.exist("system"))
    reconstruction.use_system_matrix();
  if (!cl.exist("system")) {
    reconstruction.calculate_weight();
  }

  reconstruction.calculate_sensitivity();

  {
    std::ofstream out_sensitivity(output.wo_ext() + "_sensitivity" +
                                  output.ext());
    for (auto& sens : reconstruction.sensitivity())
      out_sensitivity << sens << "\n";
  }

  std::ifstream response_stream(cl.get<std::string>("response"));
  reconstruction.fscanf_responses(response_stream);
  if (verbose)
    std::cout << "read in responses\n";

  if (cl.exist("graphics")) {
    int event_num = cl.get<int>("event");
    auto graphics_file_name = output.wo_ext() + ".m";
    std::ofstream out_graphics(graphics_file_name);

    MathematicaGraphics graphics(out_graphics);

    auto scanner = ScannerBuilder::build_multiple_rings(
        PET2D_BARREL_SCANNER_CL(cl, SquareDetector::F));
    graphics.add(scanner);

    auto event = reconstruction.event(event_num);
    auto lor = event.lor;
    graphics.add(scanner, lor);
    graphics.add(event.p);
    for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
      graphics.add_pixel(geometry.grid, it->pixel);
    }

    return 0;
  }

  auto n_blocks = cl.get<int>("blocks");
  auto n_iterations = cl.get<int>("iterations");

  const int n_file_digits = n_blocks * n_iterations >= 1000
                                ? 4
                                : n_blocks * n_iterations >= 100
                                      ? 3
                                      : n_blocks * n_iterations >= 10 ? 2 : 1;

  util::progress progress(verbose, n_blocks * n_iterations, 1);

  for (int block = 0; block < n_blocks; ++block) {
    for (int i = 0; i < n_iterations; i++) {
      progress(block * n_iterations + i);
      reconstruction();
      progress(block * n_iterations + i, true);
    }
    std::stringstream fn;
    fn << output_base_name << "_"      // phantom_
       << std::setw(n_file_digits)     //
       << std::setfill('0')            //
       << (block + 1) * n_iterations;  // 001
    util::obstream out(fn.str() + ".bin");
    out << reconstruction.rho();
  }

  CMDLINE_CATCH
}

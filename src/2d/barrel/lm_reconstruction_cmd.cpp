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
#include "2d/barrel/lor_geometry.h"
#include "2d/barrel/sparse_matrix.h"
#include "2d/barrel/lm_reconstruction.h"

#include "common/types.h"
#include "common/mathematica_graphics.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/lm_reconstruction.h"
#include "simple_geometry.h"
#endif

using RNG = std::mt19937;
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<Scanner>;
using MathematicaGraphics = Common::MathematicaGraphics<F>;
#if HAVE_CUDA
using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
#endif

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET2D::Barrel::add_lm_reconstruction_options(cl);
  cl.add("events", '\0', "print events");
  cl.parse_check(argc, argv);
  PET2D::Barrel::calculate_scanner_options(cl, argc);

  auto verbose = cl.count("verbose");

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  util::ibstream in_geometry(cl.get<cmdline::path>("geometry"));
  if (!in_geometry.is_open()) {
    throw("cannot open geometry file: " + cl.get<cmdline::path>("geometry"));
  }
  PET2D::Barrel::Geometry<F, S> geometry(in_geometry);

  if (verbose) {
    std::cout << "LM reconstruction:" << std::endl
              << "   detectors = " << geometry.n_detectors << std::endl
              << "  pixel_grid = " << geometry.grid.n_columns << " x "
              << geometry.grid.n_rows << " / " << geometry.grid.pixel_size
              << std::endl;
  }

  if (cl.exist("system")) {
    geometry.erase_pixel_info();
    auto system_matrix_file_name = cl.get<cmdline::path>("system");
    util::ibstream in_matrix(system_matrix_file_name);
    if (!in_matrix.is_open()) {
      throw("cannot open system matrix file: " +
            cl.get<cmdline::path>("system"));
    }
    PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit> matrix(in_matrix);
    if (verbose) {
      std::cout << "read in system matrix: " << cl.get<cmdline::path>("system")
                << std::endl;
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

  if (cl.exist("system")) {
    reconstruction.use_system_matrix();
  } else {
    reconstruction.calculate_weight();
  }
  reconstruction.calculate_sensitivity();

  for (const auto& fn : cl.rest()) {
    std::ifstream in_response(fn);
    if (!in_response.is_open()) {
      throw("cannot open response file: " + fn);
    }
    reconstruction << in_response;
  }
  if (verbose) {
    std::cout << "      events = " << reconstruction.n_events() << std::endl;
  }

  auto output_name = cl.get<cmdline::path>("output");
  auto output_base_name = output_name.wo_ext();
  auto output_ext = output_name.ext();
  auto output_txt = output_ext == ".txt";

  if (cl.exist("graphics")) {
    int event_num = cl.get<int>("event");
    auto graphics_file_name = output_base_name + ".m";
    std::ofstream out_graphics(graphics_file_name);

    MathematicaGraphics graphics(out_graphics);

    auto scanner = ScannerBuilder::build_multiple_rings(
        PET2D_BARREL_SCANNER_CL(cl, SquareDetector::F));
    graphics.add(scanner);

    auto event = reconstruction.event(event_num);
    auto lor = event.lor;
    graphics.add(scanner, lor);
#if FULL_EVENT_INFO
    graphics.add(event.p);
#endif
    const auto& lor_geometry = geometry[event.lor];
    for (auto i = event.first_pixel_info_index; i < event.last_pixel_info_index;
         ++i) {
      const auto& pixel_info = lor_geometry.pixel_infos[i];
      graphics.add_pixel(geometry.grid, pixel_info.pixel);
    }

    return 0;
  }

  if (cl.exist("events")) {
    std::ofstream out(output_base_name + "_events.txt");
    for (size_t i = 0; i < reconstruction.n_events(); i++) {
      out << reconstruction.event(i).p << " " << reconstruction.event(i).t
          << "\n";
    }
    return 0;
  }

  auto n_blocks = cl.get<int>("blocks");
  auto n_iterations_in_block = cl.get<int>("iterations");
  auto n_iterations = n_blocks * n_iterations_in_block;

  util::progress progress(verbose, n_iterations, 1);

#if HAVE_CUDA
  if (cl.exist("gpu")) {
    SimpleGeometry simple_geometry(geometry);
    PET2D::Barrel::GPU::LMReconstruction::run(
        simple_geometry,
        reconstruction.events().data(),
        reconstruction.n_events(),
        reconstruction.sigma(),
        geometry.grid.n_columns,
        geometry.grid.n_rows,
        n_blocks,
        n_iterations_in_block,
        [&](int iteration,
            const PET2D::Barrel::GPU::LMReconstruction::Output& output) {
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
    if (output_base_name.length()) {
      std::ofstream out_sensitivity(output_base_name + "_sensitivity" +
                                    output_ext);
      out_sensitivity << reconstruction.sensitivity();
      util::png_writer png_sensitivity(output_base_name + "_sensitivity.png");
      png_sensitivity << reconstruction.sensitivity();
    }

    for (int block = 0; block < n_blocks; ++block) {
      for (int i = 0; i < n_iterations_in_block; i++) {
        progress(block * n_iterations_in_block + i);
        reconstruction();
        progress(block * n_iterations_in_block + i, true);
      }
      if (!output_base_name.length())
        continue;
      auto fn = output_base_name.add_index((block + 1) * n_iterations_in_block,
                                           n_iterations);
      util::png_writer png(fn + ".png");
      png << reconstruction.rho();
      if (output_txt) {
        std::ofstream txt(fn + ".txt");
        txt << reconstruction.rho();
      } else {
        util::obstream bin(fn + output_ext);
        bin << reconstruction.rho();
      }
    }
  }

  CMDLINE_CATCH
}

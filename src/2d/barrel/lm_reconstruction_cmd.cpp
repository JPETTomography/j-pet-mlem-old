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
using Output = PET2D::Barrel::LMReconstruction<F, S>::Output;
using Geometry = PET2D::Barrel::Geometry<F, S>;
using Reconstruction = PET2D::Barrel::LMReconstruction<F, S>;
#if HAVE_CUDA
using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
#endif

void print_geometry_info(const Geometry& geometry) {
  std::cout << "LM reconstruction:" << std::endl
            << "   detectors = " << geometry.n_detectors << std::endl
            << "  pixel_grid = " << geometry.grid.n_columns << " x "
            << geometry.grid.n_rows << " / " << geometry.grid.pixel_size
            << std::endl;
}

void read_system_matrix_into_geometry(std::string system_matrix_file_name,
                                      Geometry& geometry,
                                      bool verbose = false) {
  geometry.erase_pixel_info();

  util::ibstream in_matrix(system_matrix_file_name);
  if (!in_matrix.is_open()) {
    throw("cannot open system matrix file: " + system_matrix_file_name);
  }
  PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit> matrix(in_matrix);
  if (verbose) {
    std::cout << "read in system matrix: " << system_matrix_file_name
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

void output_sensitivity(const cmdline::path& output,
                        const Output& sensitivity) {
  std::ofstream out_sensitivity(output.wo_ext() + "_sensitivity" +
                                output.ext());
  out_sensitivity << sensitivity;
  util::png_writer png_sensitivity(output.wo_ext() + "_sensitivity.png");
  png_sensitivity << sensitivity;
}

void read_in_response_files(const cmdline::parser& cl,
                            Reconstruction& reconstruction,
                            bool verbose = false) {
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
}

void plot_event_to_mathematica(int event_num,
                               const cmdline::parser cl,
                               const cmdline::path& output,
                               const Geometry& geometry,
                               const Reconstruction& reconstruction) {

  auto graphics_file_name = output.wo_ext() + ".m";
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
}

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET2D::Barrel::add_lm_reconstruction_options(cl);
  cl.parse_check(argc, argv);
  PET2D::Barrel::calculate_scanner_options(cl, argc);

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
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
    print_geometry_info(geometry);
  }

  if (cl.exist("system")) {
    auto system_matrix_file_name = cl.get<cmdline::path>("system");
    read_system_matrix_into_geometry(
        system_matrix_file_name, geometry, verbose);
  }

  PET2D::Barrel::LMReconstruction<F, S> reconstruction(
      geometry, cl.get<double>("s-dl") / 2);

  if (cl.exist("system")) {
    reconstruction.use_system_matrix();
  } else {
    reconstruction.calculate_weight();
  }

  reconstruction.calculate_sensitivity();

  output_sensitivity(output, reconstruction.sensitivity());

  read_in_response_files(cl, reconstruction, verbose);

  if (cl.exist("graphics")) {
    int event_num = cl.get<int>("event");
    plot_event_to_mathematica(event_num, cl, output, geometry, reconstruction);
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
        [&](int iteration, float* output) {
          auto fn = iteration >= 0
                        ? output_base_name.add_index(iteration, n_iterations)
                        : output_base_name + "_sensitivity";
          util::png_writer png(fn + ".png");
          png.write(geometry.grid.n_columns, geometry.grid.n_rows, output);
          std::ofstream txt(fn + ".txt");
          for (int i = 0; i < geometry.grid.n_columns * geometry.grid.n_rows;
               ++i) {
            txt << output[i];
            if ((i + 1) % geometry.grid.n_columns == 0) {
              txt << "\n";
            } else {
              txt << " ";
            }
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
        progress(block * n_iterations_in_block + i);
        reconstruction();
        progress(block * n_iterations_in_block + i, true);
      }
      auto fn = output_base_name.add_index((block + 1) * n_iterations_in_block,
                                           n_iterations);
      util::obstream out(fn + ".bin");
      out << reconstruction.rho().as_vector();
      util::png_writer png(fn + ".png");
      png << reconstruction.rho();
    }
  }

  CMDLINE_CATCH
}

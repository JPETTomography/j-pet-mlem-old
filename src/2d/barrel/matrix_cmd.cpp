/// \page cmd_2d_barrel_matrix 2d_barrel_matrix
/// \brief 2D Barrel PET system matrix construction tool
///
/// Creates system matrix file and accomanying SVG & PNG files for
/// reconstruction \ref cmd_2d_barrel_reconstruction.
///
/// Example
/// -------
///
/// - Create system matrix for 2 rings of 48 detectors using 1 million emissions
///   from each pixel:
///
///        2d_barrel_matrix -s square -w 0.007 -h 0.017
///                         -r 0.360 -d 48 --radius2 0.400
///                         -e 1000000 -o data/201412_rings/gpu_2rings
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/matrix_cmd.txt
///
/// \sa \ref cmd_2d_barrel_phantom, \ref cmd_2d_barrel_reconstruction

#include <iostream>
#include <random>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/random.h"
#include "util/png_writer.h"
#include "util/svg_ostream.h"
#include "util/progress.h"
#include "util/variant.h"
#include "util/backtrace.h"

#include "scanner_builder.h"
#include "ring_scanner.h"
#include "generic_scanner.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"
#include "matrix_pixel_major.h"

#include "2d/geometry/pixel.h"
#include "lor.h"
#include "common/model.h"

#include "options.h"
#include "monte_carlo.h"

#include "common/types.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/matrix.h"
#endif

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;

template <class DetectorClass>
using Scanner = PET2D::Barrel::GenericScanner<DetectorClass, S>;

// all available detector shapes
using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;
using CircleScanner = Scanner<PET2D::Barrel::CircleDetector<F>>;
using TriangleScanner = Scanner<PET2D::Barrel::TriangleDetector<F>>;
using HexagonalScanner = Scanner<PET2D::Barrel::PolygonalDetector<6, F>>;

template <class ScannerClass, class ModelClass>
void print_parameters(cmdline::parser& cl, const ScannerClass& scanner);

using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
using ComputeMatrix = PET2D::Barrel::MatrixPixelMajor<Pixel, LOR, Hit>;

template <class DetectorClass, class ModelClass>
static SparseMatrix run(cmdline::parser& cl,
                        DetectorClass& scanner,
                        ModelClass& model);

template <class ScannerClass>
void post_process(cmdline::parser& cl,
                  ScannerClass& scanner,
                  SparseMatrix& sparse_matrix);

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET2D::Barrel::add_matrix_options(cl);
  cl.parse_check(argc, argv);
  cmdline::load_accompanying_config(cl, false);
  PET2D::Barrel::calculate_scanner_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  if (cl.exist("png") && !cl.exist("from")) {
    throw("need to specify --from lor when output --png option is specified");
  }
  if (!cl.exist("png") && cl.exist("from")) {
    throw("need to specify output --png option when --from is specified");
  }

  const auto& shape = cl.get<std::string>("shape");
  const auto& model_name = cl.get<std::string>("model");
  const auto& length_scale = cl.get<double>("base-length");

// these are wrappers running actual simulation
#if HAVE_CUDA
#define _RUN(cl, scanner, model)                        \
  cl.exist("gpu") ? PET2D::Barrel::GPU::Matrix::run(cl) \
                  : run(cl, scanner, model)
#else
#define _RUN(cl, scanner, model) run(cl, scanner, model)
#endif
#define RUN(detector_type, model_type, ...)                               \
  detector_type scanner =                                                 \
      PET2D::Barrel::ScannerBuilder<detector_type>::build_multiple_rings( \
          PET2D_BARREL_SCANNER_CL(cl, detector_type::F));                 \
  model_type model __VA_ARGS__;                                           \
  print_parameters<detector_type, model_type>(cl, scanner);               \
  auto sparse_matrix = _RUN(cl, scanner, model);                          \
  post_process(cl, scanner, sparse_matrix)

  // run simmulation on given detector model & shape
  if (model_name == "always") {
    if (shape == "square") {
      RUN(SquareScanner, Common::AlwaysAccept<F>);
    } else if (shape == "circle") {
      RUN(CircleScanner, Common::AlwaysAccept<F>);
    } else if (shape == "triangle") {
      RUN(TriangleScanner, Common::AlwaysAccept<F>);
    } else if (shape == "hexagon") {
      RUN(HexagonalScanner, Common::AlwaysAccept<F>);
    }
  } else if (model_name == "scintillator") {
    if (shape == "square") {
      RUN(SquareScanner, Common::ScintillatorAccept<F>, (length_scale));
    } else if (shape == "circle") {
      RUN(CircleScanner, Common::ScintillatorAccept<F>, (length_scale));
    } else if (shape == "triangle") {
      RUN(TriangleScanner, Common::ScintillatorAccept<F>, (length_scale));
    } else if (shape == "hexagon") {
      RUN(HexagonalScanner, Common::ScintillatorAccept<F>, (length_scale));
    }
  }

  CMDLINE_CATCH
}

template <class ScannerClass, class ModelClass>
void print_parameters(cmdline::parser& cl, const ScannerClass& scanner) {
  auto& n_pixels = cl.get<int>("n-pixels");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.count("verbose");
  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = 0;
    n_tof_positions = scanner.n_tof_positions(tof_step, max_bias);
  }
  if (verbose) {
    std::cerr << "Monte-Carlo:" << std::endl;
#if _OPENMP
    std::cerr << "   OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif
    std::cerr << "    pixels in row = " << n_pixels << std::endl;
    std::cerr << "       pixel size = " << s_pixel << std::endl;
    std::cerr << "       fov radius = " << scanner.fov_radius() << std::endl;
    std::cerr << "     outer radius = " << scanner.outer_radius() << std::endl;
    std::cerr << "         max bias = " << max_bias << std::endl;
    std::cerr << "         TOF step = " << tof_step << std::endl;
    std::cerr << "    TOF positions = " << n_tof_positions << std::endl;
    std::cerr << "        emissions = " << n_emissions << std::endl;
  }
}

template <class DetectorClass, class ModelClass>
static SparseMatrix run(cmdline::parser& cl,
                        DetectorClass& scanner,
                        ModelClass& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.count("verbose");

  std::random_device rd;
  util::random::tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = 0;
    n_tof_positions = scanner.n_tof_positions(tof_step, max_bias);
  }

  ComputeMatrix::SparseMatrix sparse_matrix(
      n_pixels, scanner.size(), n_tof_positions);

  for (auto& fn : cl.rest()) {
    util::ibstream in(fn, std::ios::binary);
    if (!in.is_open())
      throw("cannot open input file: " + fn);
    try {
      ComputeMatrix::SparseMatrix in_sparse_matrix(in);
      if (verbose) {
        std::cerr << "read sparse matrix: " << fn << std::endl;
        std::cerr << " pixels in row = " << in_sparse_matrix.n_pixels_in_row()
                  << std::endl;
        std::cerr << " pixel size     = " << cl.get<double>("s-pixel");
        std::cerr << " TOF positions = " << in_sparse_matrix.n_tof_positions()
                  << std::endl;
        std::cerr << " emissions     = " << in_sparse_matrix.n_emissions()
                  << std::endl;
        std::cerr << std::endl;
      }
      if (sparse_matrix.empty()) {
        sparse_matrix = in_sparse_matrix;
        // if we don't have stuff set, set it using matrix
        if (!cl.exist("n-pixels"))
          n_pixels = sparse_matrix.n_pixels_in_row();
        if (!cl.exist("tof-step")) {
          n_tof_positions = sparse_matrix.n_tof_positions();
          if (n_emissions > 0) {
            throw("TOF step must be specified for input TOF matrix: " + fn);
          }
        }
      } else {
        // join with previous matrix
        sparse_matrix << in_sparse_matrix;
      }
    } catch (std::string& ex) {
      throw(ex + ": " + fn);
    } catch (const char* ex) {
      throw(std::string(ex) + ": " + fn);
    }
  }

  ComputeMatrix matrix(n_pixels, scanner.size(), n_tof_positions);
  if (!sparse_matrix.empty()) {
    matrix << sparse_matrix;
    sparse_matrix.resize(0);
  }

#ifdef __linux__
  struct timespec start, stop;
  clock_gettime(CLOCK_REALTIME, &start);
#endif

  PET2D::Barrel::MonteCarlo<DetectorClass, ComputeMatrix> monte_carlo(
      scanner, matrix, s_pixel, tof_step, m_pixel);
  util::progress progress(verbose, matrix.total_n_pixels_in_triangle, 1);
  monte_carlo(gen, model, n_emissions, progress);

#ifdef GPU_TOF_TEST
  monte_carlo.test(gen, model, n_emissions);
#endif

#ifdef __linux__
  if (verbose) {
    clock_gettime(CLOCK_REALTIME, &stop);
    std::cerr << "time : "
              << ((1.0e9 * stop.tv_sec + stop.tv_nsec) -
                  (1.0e9 * start.tv_sec + start.tv_nsec)) /
                     1.0e9
              << std::endl;
  }
#endif

  return matrix.to_sparse();
}

template <class ScannerClass>
void post_process(cmdline::parser& cl,
                  ScannerClass& scanner,
                  SparseMatrix& sparse_matrix) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_detectors = cl.get<int>("n-detectors");

  // generate output
  if (cl.exist("output")) {
    auto fn = cl.get<cmdline::path>("output");
    auto fn_wo_ext = fn.wo_ext();
    auto fn_wo_path = fn_wo_ext.wo_path();
    bool full = cl.exist("full");
    util::obstream out(fn, std::ios::binary | std::ios::trunc);
    if (full) {
      auto full_matrix = sparse_matrix.to_full(scanner.symmetry_descriptor());
      out << full_matrix;
    } else {
      out << sparse_matrix;
    }

    std::ofstream os(fn_wo_ext + ".cfg", std::ios::trunc);
    os << cl;

    try {
      util::png_writer png(fn_wo_ext + ".png");
      sparse_matrix.output_bitmap(png);
    } catch (const char* ex) {
      // don't bail out just produce warning
      std::cerr << "warning: " << ex << std::endl;
    }

    util::svg_ostream<F> svg(fn_wo_ext + ".svg",
                             scanner.outer_radius(),
                             scanner.outer_radius(),
                             1024.,
                             1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << scanner;
  }

  // visual debugging output
  if (cl.exist("png")) {
    LOR lor(0, 0);
    lor.first = cl.get<int>("from");
    if (cl.exist("to")) {
      lor.second = cl.get<int>("to");
    } else {
      lor.second = (lor.first + n_detectors / 2) % n_detectors;
    }
    if (lor.first < lor.second)
      std::swap(lor.first, lor.second);

    auto fn = cl.get<cmdline::path>("png");
    auto fn_wo_ext = fn.wo_ext();
    auto fn_wo_path = fn_wo_ext.wo_path();

    util::png_writer png(fn);
    auto position = cl.get<int>("pos");
    if (cl.exist("full")) {
      sparse_matrix.to_full(scanner.symmetry_descriptor())
          .output_bitmap(png, lor, position);
    } else {
      sparse_matrix.output_bitmap(png, lor, position);
    }

    util::svg_ostream<F> svg(fn_wo_ext + ".svg",
                             scanner.outer_radius(),
                             scanner.outer_radius(),
                             1024.,
                             1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << scanner;
  }
}

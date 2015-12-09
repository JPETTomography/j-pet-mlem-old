/// \page cmd_3d_hybrid_matrix 3d_hybrid_matrix
/// \brief 3D Longitudinal PET system matrix construction tool
///
/// Creates system matrix file and accomanying SVG & PNG files for single slice
/// along z-axis.
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verboutput 3d_hybrid_matrix
///
/// \sa \ref cmd_3d_hybrid_phantom

#include <random>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/sparse_matrix.h"
#include "2d/barrel/matrix_pixel_major.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "common/model.h"
#include "monte_carlo.h"
#include "util/progress.h"

#include "util/png_writer.h"
#include "util/svg_ostream.h"
#include "util/json.h"
#include "util/backtrace.h"

#include "options.h"

#include "common/types.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/matrix.h"
#endif

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;

using Pixel = PET2D::Pixel<Scanner2D::S>;
using LOR = Scanner2D::LOR;
using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
using ComputeMatrix = PET2D::Barrel::MatrixPixelMajor<Pixel, LOR, Hit>;

template <class ScannerClass, class ModelClass, typename... ModelArgs>
static void run(cmdline::parser& cl, ModelArgs... args);

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Hybrid::add_matrix_options(cl);
  cl.parse_check(argc, argv);
  cmdline::load_accompanying_config(cl, false);
  if (!cl.exist("tof-step")) {
    cl.dontsave("tof-step"), cl.dontsave("s-dl");
  }
  PET3D::Hybrid::calculate_scanner_options(cl, argc);

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
  if (shape != "square") {
    throw("only square is supported, unsupported shape: " + shape);
  }

  const auto& model_name = cl.get<std::string>("model");
  const auto& length_scale = cl.get<double>("base-length");

  if (model_name == "always") {
    run<Scanner2D, Common::AlwaysAccept<F>>(cl);
  } else if (model_name == "scintillator") {
    run<Scanner2D, Common::ScintillatorAccept<F>>(cl, length_scale);
  } else {
    throw("unknown model: " + model_name);
  }

  CMDLINE_CATCH
}

template <class ScannerClass, class ModelClass, typename... ModelArgs>
static void run(cmdline::parser& cl, ModelArgs... args) {

  auto scanner2D =
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, Scanner::F));
  Scanner scanner(scanner2D, F(cl.get<double>("length")));
  ModelClass model(args...);

  auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto verbose = cl.count("verbose");
  auto& z_position = cl.get<double>("z-position");

  if (verbose) {
    std::cerr << "  n pixels = " << n_pixels << std::endl;
    std::cerr << "pixel size = " << s_pixel << std::endl;
  }

  std::random_device rd;
  util::random::tausworthe rng(rd());
  if (cl.exist("seed")) {
    rng.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  if (verbose) {
    std::cerr << "Monte-Carlo:" << std::endl;
#if _OPENMP
    if (!cl.exist("gpu"))
      std::cerr << " OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif
    std::cerr << "  pixels in row = " << n_pixels << std::endl;
    std::cerr << "     pixel size = " << s_pixel << std::endl;
    std::cerr << "     fov radius = " << scanner2D.fov_radius() << std::endl;
    std::cerr << "   outer radius = " << scanner2D.outer_radius() << std::endl;
    std::cerr << "      emissions = " << n_emissions << std::endl;
  }

  ComputeMatrix::SparseMatrix sparse_matrix(n_pixels, scanner.barrel.size());

  for (auto& fn : cl.rest()) {
    util::ibstream in(fn, std::ios::binary);
    if (!in.is_open())
      throw("cannot open input file: " + fn);
    try {
      ComputeMatrix::SparseMatrix in_sparse_matrix(in);
      if (in_sparse_matrix.n_tof_positions() > 1)
        throw("hybrid Monte-Carlo does not support TOF positions");
      if (verbose) {
        std::cerr << "read sparse matrix: " << fn << std::endl;
        std::cerr << "  pixels in row = " << in_sparse_matrix.n_pixels_in_row()
                  << std::endl;
        std::cerr << "     pixel size = " << cl.get<double>("s-pixel");
        std::cerr << "      emissions = " << in_sparse_matrix.n_emissions()
                  << std::endl;
        std::cerr << std::endl;
      }
      if (sparse_matrix.empty()) {
        sparse_matrix = in_sparse_matrix;
        // if we don't have stuff set, set it using matrix
        if (!cl.exist("n-pixels"))
          n_pixels = sparse_matrix.n_pixels_in_row();
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

  ComputeMatrix matrix(n_pixels, scanner.barrel.size());
  util::progress progress(verbose, matrix.total_n_pixels_in_triangle, 1);

  if (!n_emissions) {
    // don't run any simulation computation
  }
#if HAVE_CUDA
  // GPU computation
  else if (cl.exist("gpu")) {
    bool was_empty = sparse_matrix.empty();
    PET3D::Hybrid::GPU::Matrix::run<Scanner>(
        scanner,
        rng,
        n_emissions,
        z_position,
        n_pixels,
        s_pixel,
        cl.get<double>("base-length"),
        [&](int completed, bool finished) { progress(completed, finished); },
        [&](LOR lor, Pixel pixel, Hit hits) {
          sparse_matrix.emplace_back(lor, 0, pixel, hits);
        },
        cl.get<int>("cuda-device"),
        cl.get<int>("cuda-blocks"),
        cl.get<int>("cuda-threads"),
        [&](const char* device_name, int n_final_emissions) {
          n_emissions = n_final_emissions;
          if (verbose) {
            std::cerr << "    CUDA device = " << device_name << std::endl;
            std::cerr << "final emissions = " << n_final_emissions << std::endl;
          }
        });
    sparse_matrix.increment_n_emissions(n_emissions);
    if (!was_empty) {
      sparse_matrix.sort_by_lor_n_pixel();
      sparse_matrix.merge_duplicates();
    }
  }
#endif
  // CPU reference computation
  else {
    if (!sparse_matrix.empty()) {
      matrix << sparse_matrix;
      sparse_matrix.resize(0);
    }
    PET3D::Hybrid::MonteCarlo<Scanner, ComputeMatrix> monte_carlo(
        scanner, matrix, s_pixel, m_pixel);
    monte_carlo(z_position, rng, model, n_emissions, progress);
    sparse_matrix = matrix.to_sparse();
  }

  // generate output
  if (cl.exist("output")) {
    auto fn = cl.get<cmdline::path>("output");
    auto fn_wo_ext = fn.wo_ext();
    auto fn_wo_path = fn_wo_ext.wo_path();
    bool full = cl.exist("full");
    util::obstream out(fn, std::ios::binary | std::ios::trunc);
    if (full) {
      auto full_matrix =
          sparse_matrix.to_full(scanner.barrel.symmetry_descriptor());
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
                             scanner.barrel.outer_radius(),
                             scanner.barrel.outer_radius(),
                             1024.,
                             1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << const_cast<Scanner2D&>(scanner.barrel);
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
      sparse_matrix.to_full(scanner.barrel.symmetry_descriptor())
          .output_bitmap(png, lor, position);
    } else {
      sparse_matrix.output_bitmap(png, lor, position);
    }

    util::svg_ostream<F> svg(fn_wo_ext + ".svg",
                             scanner.barrel.outer_radius(),
                             scanner.barrel.outer_radius(),
                             1024.,
                             1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << const_cast<Scanner2D&>(scanner.barrel);
  }
}

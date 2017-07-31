/// \page cmd_2d_barrel_matrix 2d_barrel_matrix
/// \brief 2D Barrel PET system matrix construction tool
///
/// Creates system matrix file and accomanying SVG & PNG files for
/// reconstruction \ref cmd_2d_barrel_reconstruction.
///
/// Example
/// -------
///
/// 1. Make a \c playground directory and step into it
///
///        mkdir playground
///        cd playground
///
/// 2. Create system matrix for 2 rings of 36cm and 40cm radius (\c -r)
///    and 48 detectors (\c -d) using 1 million emissions
///    from each pixel (\c -e) of size 4mm (\c -p) and be verbose (\c -v)
///
///        ../2d_barrel_matrix \
///          -s square -w 0.007 -h 0.017 \
///          -r 0.360,0.400 -d 48 \
///          -e 1m \
///          -p 0.004 \
///          -o m_2rings \
///          -v
///
/// 3. Create full system matrix for "big" barrel configuration
///    1 million emissions from each pixel (\c -e) of size 4mm (\c -p)
///    and be verbose (\c -v)
///
///        ../2d_barrel_matrix \
///          -c ../config/big.cfg \
///          -e 1m \
///          -p 0.004 \
///          --full \
///          -o f_big \
///          -v
///
/// \note Unless \c --full is given, the command produces triangular (1/8)
/// matrix which is unsuitable for reconstruction.
///
/// \note Generating 1 million emissions system matrix may take a while on CPU,
/// it is recommended to use a CUDA compatible (NVIDIA) GPU card enabled by
/// passing \c -G or \c --gpu option.
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verboutput 2d_barrel_matrix
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

#include "2d/barrel/symmetry_descriptor.h"

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

template <class ScannerClass>
using Scanner = PET2D::Barrel::GenericScanner<ScannerClass, S>;

using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;
using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
using ComputeMatrix = PET2D::Barrel::MatrixPixelMajor<Pixel, LOR, Hit>;

template <class ScannerClass, class ModelClass, typename... ModelArgs>
static void run(cmdline::parser& cl, ModelArgs... args);

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;

  std::ostringstream msg;
  msg << "matrix_file ..." << std::endl;
  msg << "build: " << VARIANT << std::endl;
  msg << "note: All length options below should be expressed in meters.";
  cl.footer(msg.str());

  PET2D::Barrel::add_config_option(cl);
  PET2D::Barrel::add_pixel_options(cl);
  // PET2D::Barrel::add_scanner_options(cl);
  PET2D::Barrel::add_model_options(cl);

  cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
  cl.add<size_t>("n-emissions",
                 'e',
                 "emissions per pixel",
                 false,
                 0,
                 cmdline::not_from_file);
  cl.add<cmdline::path>("output",
                        'o',
                        "output binary triangular/full sparse system matrix",
                        cmdline::dontsave);
  cl.add("full", 'f', "output full non-triangular sparse system matrix");

  // visual debugging params
  cl.add<cmdline::path>("png", 0, "output lor to png", cmdline::dontsave);
  cl.add<int>("from", 0, "lor start detector to output", cmdline::dontsave, -1);
  cl.add<int>("to", 0, "lor end detector to output", cmdline::dontsave, -1);
  cl.add<int>("pos", 0, "position to output", cmdline::dontsave, -1);
  // printing & stats params
  cl.add("print", 0, "print triangular sparse system matrix");
  cl.add("stats", 0, "show stats");
  cl.add("wait", 0, "wait before exit");
  cl.add("verbose", 'v', "prints the iterations information on std::out");
  cl.add<util::random::tausworthe::seed_type>(
      "seed", 'S', "random number generator seed", cmdline::dontsave);
  Common::add_cuda_options(cl);
  Common::add_openmp_options(cl);

  cl.add<cmdline::path>(
      "detector-file", '\0', "detector description file", false);
  cl.add<cmdline::path>(
      "detector-file-sym", '\0', "detector symmetries description file", false);

  cl.parse_check(argc, argv);
  cmdline::load_accompanying_config(cl, false);

  if (!cl.exist("tof-step")) {
    cl.dontsave("tof-step"), cl.dontsave("s-dl");
  }

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

  const auto& model_name = cl.get<std::string>("model");
  const auto& length_scale = cl.get<double>("base-length");

  // run simmulation on given detector model & shape
  if (model_name == "always")
    run<SquareScanner, Common::AlwaysAccept<F>>(cl);
  else if (model_name == "scintillator")
    run<SquareScanner, Common::ScintillatorAccept<F>>(cl, length_scale);
  else {
    std::cerr << "model `" << model_name << "' is not supported\n";
    exit(1);
  }

  CMDLINE_CATCH
}

template <class ScannerClass, class ModelClass, typename... ModelArgs>
static void run(cmdline::parser& cl, ModelArgs... args) {

  std::ifstream in_dets(cl.get<cmdline::path>("detector-file"));
  auto scanner =
      PET2D::Barrel::ScannerBuilder<ScannerClass>::deserialize(in_dets);
  in_dets.close();

  std::ifstream in_syms(cl.get<cmdline::path>("detector-file-sym"));
  auto symmetry = PET2D::Barrel::SymmetryDescriptor<S>::deserialize(in_syms);
  scanner.set_symmetry_descriptor(symmetry);
  in_syms.close();

  ModelClass model(args...);

  auto n_detectors = scanner.size();
  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<size_t>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.count("verbose");

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = 0;
    n_tof_positions = scanner.n_tof_positions(tof_step, max_bias);
  }

  std::random_device rd;
  util::random::tausworthe rng(rd());
  if (cl.exist("seed")) {
    rng.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  if (verbose) {
    std::cerr << "Monte-Carlo:" << std::endl;
#if _OPENMP
#if HAVE_CUDA
    if (!cl.exist("gpu"))
#endif
      std::cerr << " OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif
    std::cerr << "  pixels in row = " << n_pixels << std::endl;
    std::cerr << "     pixel size = " << s_pixel << std::endl;
    std::cerr << "     fov radius = " << scanner.fov_radius() << std::endl;
    std::cerr << "   outer radius = " << scanner.outer_radius() << std::endl;
    std::cerr << "       max bias = " << max_bias << std::endl;
    std::cerr << "       TOF step = " << tof_step << std::endl;
    std::cerr << "  TOF positions = " << n_tof_positions << std::endl;
    std::cerr << "      emissions = " << n_emissions << std::endl;
  }

  ComputeMatrix::SparseMatrix sparse_matrix(
      n_pixels, scanner.size(), n_tof_positions);

  for (auto& fn : cl.rest()) {
    util::ibstream in(fn, std::ios::binary);
    ENSURE_IS_OPEN(in, "input matrix", fn);
    try {
      ComputeMatrix::SparseMatrix in_sparse_matrix(in);
      if (verbose) {
        std::cerr << "read sparse matrix: " << fn << std::endl;
        std::cerr << "  pixels in row = " << in_sparse_matrix.n_pixels_in_row()
                  << std::endl;
        std::cerr << "     pixel size = " << cl.get<double>("s-pixel");
        std::cerr << "  TOF positions = " << in_sparse_matrix.n_tof_positions()
                  << std::endl;
        std::cerr << "      emissions = " << in_sparse_matrix.n_emissions()
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
  util::progress progress(verbose, matrix.total_n_pixels_in_triangle, 1);

  if (!n_emissions) {
    // don't run any simulation computation
  }
#if HAVE_CUDA
  // GPU computation
  else if (cl.exist("gpu")) {
    bool was_empty = sparse_matrix.empty();
    PET2D::Barrel::GPU::Matrix::run<ScannerClass>(
        scanner,
        rng,
        n_emissions,
        tof_step,
        n_tof_positions,
        n_pixels,
        s_pixel,
        cl.get<double>("base-length"),
        [&](int completed, bool finished) { progress(completed, finished); },
        [&](LOR lor, S position, Pixel pixel, Hit hits) {
          sparse_matrix.emplace_back(lor, position, pixel, hits);
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
    PET2D::Barrel::MonteCarlo<ScannerClass, ComputeMatrix> monte_carlo(
        scanner, matrix, s_pixel, tof_step, m_pixel);
    monte_carlo(rng, model, n_emissions, progress);
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
      // lor.second = (lor.first + n_detectors[0] / 2) % n_detectors[0];
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

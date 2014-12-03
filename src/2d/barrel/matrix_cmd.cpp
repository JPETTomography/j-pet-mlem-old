/// \page _2d_barrel_matrix 2d_barrel_matrix
/// \brief 2D Barrel PET system matrix construction tool
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/matrix_cmd.txt

#include <iostream>
#include <random>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "util/random.h"
#include "detector_ring.h"
#include "detector_set.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"
#include "matrix_pixel_major.h"
#include "2d/geometry/pixel.h"
#include "lor.h"
#include "model.h"
#include "util/png_writer.h"
#include "util/svg_ostream.h"
#include "util/progress.h"
#include "util/variant.h"

#include "monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/matrix.h"
#endif

using namespace PET2D;
using namespace PET2D::Barrel;

// all available detector shapes
template <typename DetectorType>
using DetectorModel = DetectorSet<DetectorType>;
// using DetectorModel = DetectorRing<DetectorType>;
using SquareDetectorRing = DetectorModel<SquareDetector<>>;
using CircleDetectorRing = DetectorModel<CircleDetector<>>;
using TriangleDetectorRing = DetectorModel<TriangleDetector<>>;
using HexagonalDetectorRing = DetectorModel<PolygonalDetector<6>>;

template <typename DetectorRing, typename Model>
void print_parameters(cmdline::parser& cl, const DetectorRing& detector_ring);

template <typename Detector, typename Model>
static SparseMatrix<Pixel<>, LOR<>> run(cmdline::parser& cl,
                                        Detector& detector_ring,
                                        Model& model);

template <typename DetectorRing>
void post_process(cmdline::parser& cl,
                  DetectorRing& detector_ring,
                  SparseMatrix<Pixel<>, LOR<>>& sparse_matrix);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    std::ostringstream msg;
    msg << "matrix_file ..." << std::endl;
    msg << "build: " << VARIANT << std::endl;
    msg << "note: All length options below should be expressed in meters.";
    cl.footer(msg.str());

    cl.add<cmdline::path>("config",
                          'c',
                          "load config file",
                          cmdline::dontsave,
                          cmdline::path(),
                          cmdline::default_reader<cmdline::path>(),
                          cmdline::load);
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 256);
    cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
    cl.add<int>("n-detectors", 'd', "number of detectors in ring", false);
    cl.add<int>("n-detectors2", 0, "number of detectors in 2nd ring", false);
    cl.add<int>("n-detectors3", 0, "number of detectors in 3rd ring", false);
    cl.add<int>("n-emissions",
                'e',
                "emissions per pixel",
                false,
                0,
                cmdline::not_from_file);
    cl.add<double>("radius", 'r', "inner detector ring radius", false);
    cl.add<double>("ring-rotation", 0, "next ring rotation (0-1)", false);
    cl.add<double>("radius2", 0, "2nd detector ring radius", false);
    cl.add<double>("radius3", 0, "3rd detector ring radius", false);
    cl.add<double>("s-pixel", 'p', "pixel size", false);
    cl.add<double>(
        "tof-step", 't', "TOF quantisation step for distance delta", false);
    cl.add<std::string>(
        "shape",
        's',
        "detector (scintillator) shape (square, circle, triangle, hexagon)",
        false,
        "square",
        cmdline::oneof<std::string>("square", "circle", "triangle", "hexagon"));
    cl.add<double>("w-detector", 'w', "detector width", false);
    cl.add<double>("h-detector", 'h', "detector height", false);
    cl.add<double>("d-detector",
                   0,
                   "inscribe detector shape into circle of given diameter",
                   false);
    cl.add<std::string>(
        "model",
        'm',
        "acceptance model (always, scintillator)",
        false,
        "scintillator",
        cmdline::oneof<std::string>("always",
                                    "scintillator",
                                    /* obsolete */ "scintilator"));
    // NOTE: this options is obsolete (use base-length instead)
    cl.add<double>("acceptance",
                   'a',
                   "acceptance probability factor",
                   cmdline::dontsave | cmdline::hidden,
                   10.);
    cl.add<double>("base-length",
                   'l',
                   "scintillator emission base length P(l)=1-e^(-1)",
                   false,
                   0.1);
    cl.add<cmdline::path>("output",
                          'o',
                          "output binary triangular/full sparse system matrix",
                          cmdline::dontsave);
    cl.add("full", 'f', "output full non-triangular sparse system matrix");

    // visual debugging params
    cl.add<cmdline::path>("png", 0, "output lor to png", cmdline::dontsave);
    cl.add<int>(
        "from", 0, "lor start detector to output", cmdline::dontsave, -1);
    cl.add<int>("to", 0, "lor end detector to output", cmdline::dontsave, -1);
    cl.add<int>("pos", 0, "position to output", cmdline::dontsave, -1);

    // printing & stats params
    cl.add("print", 0, "print triangular sparse system matrix");
    cl.add("stats", 0, "show stats");
    cl.add("wait", 0, "wait before exit");
    cl.add("verbose", 'v', "prints the iterations information on std::out");
    cl.add<util::random::tausworthe::seed_type>(
        "seed", 'S', "random number generator seed", cmdline::dontsave);
#if HAVE_CUDA
    cl.add("gpu", 'G', "run on GPU (via CUDA)");
    cl.add<int>("cuda-blocks", 'B', "CUDA blocks", cmdline::dontsave, 64);
    cl.add<int>(
        "cuda-threads", 'W', "CUDA threads per block", cmdline::dontsave, 512);
#endif
#if _OPENMP
    cl.add<int>(
        "n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#endif

    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    // convert obsolete acceptance option to length scale
    auto& length_scale = cl.get<double>("base-length");
    if (cl.exist("acceptance") && !cl.exist("base-length")) {
      length_scale = 1.0 / cl.get<double>("acceptance");
    }
    // FIXME: fixup for spelling mistake, present in previous versions
    auto& model_name = cl.get<std::string>("model");
    if (model_name == "scintilator") {
      model_name = "scintillator";
    }

    auto& n_pixels = cl.get<int>("n-pixels");
    auto& n_detectors = cl.get<int>("n-detectors");
    auto& n_detectors2 = cl.get<int>("n-detectors2");
    auto& n_detectors3 = cl.get<int>("n-detectors3");
    auto& radius = cl.get<double>("radius");
    auto& ring_rotation = cl.get<double>("ring-rotation");
    auto& radius2 = cl.get<double>("radius2");
    auto& radius3 = cl.get<double>("radius3");
    auto& s_pixel = cl.get<double>("s-pixel");
    auto& w_detector = cl.get<double>("w-detector");
    auto& h_detector = cl.get<double>("h-detector");
    auto& d_detector = cl.get<double>("d-detector");
    auto& shape = cl.get<std::string>("shape");
    auto verbose = cl.exist("verbose");

    // check options
    if (!cl.exist("w-detector") && !cl.exist("d-detector") &&
        !cl.exist("n-detectors")) {
      throw(
          "need to specify either --w-detector, --d-detector or --n-detectors");
    }
    if (cl.exist("png") && !cl.exist("from")) {
      throw("need to specify --from lor when output --png option is specified");
    }
    if (!cl.exist("png") && cl.exist("from")) {
      throw("need to specify output --png option when --from is specified");
    }

    cmdline::load_accompanying_config(cl, false);

    if (verbose) {
      std::cerr << "assumed:" << std::endl;
    }

    // automatic radius size
    if (!cl.exist("radius")) {
      if (!cl.exist("s-pixel")) {
        radius = M_SQRT2;  // exact result
      } else {
        radius = s_pixel * n_pixels / M_SQRT2;
      }
      std::cerr << "--radius=" << radius << std::endl;
    }

    // automatic pixel size
    if (!cl.exist("s-pixel")) {
      if (!cl.exist("radius")) {
        s_pixel = 2. / n_pixels;  // exact result
      } else {
        s_pixel = M_SQRT2 * radius / n_pixels;
      }
      std::cerr << "--s-pixel=" << s_pixel << std::endl;
    }

    if (!cl.exist("w-detector")) {
      if (cl.exist("n-detectors")) {
        w_detector = 2 * M_PI * .9 * radius / n_detectors;
      } else if (cl.exist("d-detector")) {
        if (shape == "circle") {
          w_detector = d_detector;
        } else {
          auto mult = 1.;
          auto sides = 0.;
          if (shape == "triangle") {
            sides = 3.;
          } else if (shape == "square") {
            sides = 4.;
          } else if (shape == "hexagon") {
            sides = 6.;
            mult = 2.;
          } else {
            throw("cannot determine detector width for given shape");
          }
          w_detector = d_detector * std::sin(M_PI / sides) * mult;
        }
      }
      std::cerr << "--w-detector=" << w_detector << std::endl;
    }

    // automatic detector size
    // NOTE: detector height will be determined per shape
    if (!cl.exist("n-detectors")) {
      if (cl.exist("d-detector")) {
        n_detectors =
            ((int)std::floor(
                 M_PI / std::atan2(d_detector, 2 * radius + d_detector / 2)) /
             4) *
            4;
      } else {
        n_detectors =
            ((int)std::floor(M_PI / std::atan2(w_detector, 2 * radius)) / 4) *
            4;
      }
      if (!n_detectors) {
        throw("detector width is too big for given detector ring radius");
      }
      std::cerr << "--n-detectors=" << n_detectors << std::endl;
    }

// these are wrappers running actual simulation
#if HAVE_CUDA
#define _RUN(cl, detector_ring, model) \
  cl.exist("gpu") ? GPU::Matrix::run(cl) : run(cl, detector_ring, model)
#else
#define _RUN(cl, detector_ring, model) run(cl, detector_ring, model)
#endif
#define RUN(detector_type, model_type, ...)                       \
  detector_type detector_ring(n_detectors,                        \
                              radius,                             \
                              w_detector,                         \
                              h_detector,                         \
                              d_detector,                         \
                              ring_rotation,                      \
                              n_detectors2,                       \
                              radius2,                            \
                              n_detectors3,                       \
                              radius3);                           \
  model_type model{ __VA_ARGS__ };                                \
  print_parameters<detector_type, model_type>(cl, detector_ring); \
  auto sparse_matrix = _RUN(cl, detector_ring, model);            \
  post_process(cl, detector_ring, sparse_matrix)

    // run simmulation on given detector model & shape
    if (model_name == "always") {
      if (shape == "square") {
        RUN(SquareDetectorRing, AlwaysAccept<>);
      } else if (shape == "circle") {
        RUN(CircleDetectorRing, AlwaysAccept<>);
      } else if (shape == "triangle") {
        RUN(TriangleDetectorRing, AlwaysAccept<>);
      } else if (shape == "hexagon") {
        RUN(HexagonalDetectorRing, AlwaysAccept<>);
      }
    } else if (model_name == "scintillator") {
      if (shape == "square") {
        RUN(SquareDetectorRing, ScintilatorAccept<>, length_scale);
      } else if (shape == "circle") {
        RUN(CircleDetectorRing, ScintilatorAccept<>, length_scale);
      } else if (shape == "triangle") {
        RUN(TriangleDetectorRing, ScintilatorAccept<>, length_scale);
      } else if (shape == "hexagon") {
        RUN(HexagonalDetectorRing, ScintilatorAccept<>, length_scale);
      }
    }

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

template <typename DetectorRing, typename Model>
void print_parameters(cmdline::parser& cl, const DetectorRing& detector_ring) {
  auto& n_pixels = cl.get<int>("n-pixels");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.exist("verbose");
  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = detector_ring.n_tof_positions(tof_step, max_bias);
  }
  if (verbose) {
    std::cerr << "Monte-Carlo:" << std::endl;
#if _OPENMP
    std::cerr << "   OpenMP threads = " << omp_get_max_threads() << std::endl;
#endif
    std::cerr << "    pixels in row = " << n_pixels << std::endl;
    std::cerr << "     outer radius = " << detector_ring.outer_radius()
              << std::endl;
    std::cerr << "         max bias = " << max_bias << std::endl;
    std::cerr << "         TOF step = " << tof_step << std::endl;
    std::cerr << "    TOF positions = " << n_tof_positions << std::endl;
    std::cerr << "        emissions = " << n_emissions << std::endl;
  }
}

template <typename Detector, typename Model>
static SparseMatrix<Pixel<>, LOR<>> run(cmdline::parser& cl,
                                        Detector& detector_ring,
                                        Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.exist("verbose");

  std::random_device rd;
  util::random::tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = detector_ring.n_tof_positions(tof_step, max_bias);
  }

  using ComputeMatrix = MatrixPixelMajor<Pixel<>, LOR<>>;
  ComputeMatrix::SparseMatrix sparse_matrix(
      n_pixels, detector_ring.size(), n_tof_positions);

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

  ComputeMatrix matrix(n_pixels, detector_ring.size(), n_tof_positions);
  if (!sparse_matrix.empty()) {
    matrix << sparse_matrix;
    sparse_matrix.resize(0);
  }

#ifdef __linux__
  struct timespec start, stop;
  clock_gettime(CLOCK_REALTIME, &start);
#endif

  MonteCarlo<Detector, ComputeMatrix> monte_carlo(
      detector_ring, matrix, s_pixel, tof_step, m_pixel);
  util::progress progress(verbose, matrix.total_n_pixels_in_triangle(), 1);
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
                     1.0e9 << std::endl;
  }
#endif

  return matrix.to_sparse();
}

template <typename DetectorRing>
void post_process(cmdline::parser& cl,
                  DetectorRing& detector_ring,
                  SparseMatrix<Pixel<>, LOR<>>& sparse_matrix) {

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
      auto full_matrix = sparse_matrix.to_full();
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

    util::svg_ostream<> svg(fn_wo_ext + ".svg",
                            detector_ring.outer_radius(),
                            detector_ring.outer_radius(),
                            1024.,
                            1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << detector_ring;
  }

  // visual debugging output
  if (cl.exist("png")) {
    LOR<> lor(0, 0);
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
      sparse_matrix.to_full().output_bitmap(png, lor, position);
    } else {
      sparse_matrix.output_bitmap(png, lor, position);
    }

    util::svg_ostream<> svg(fn_wo_ext + ".svg",
                            detector_ring.outer_radius(),
                            detector_ring.outer_radius(),
                            1024.,
                            1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << detector_ring;
  }
}

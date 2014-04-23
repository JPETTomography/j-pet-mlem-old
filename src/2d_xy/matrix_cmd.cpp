// PET System Matrix Calculator
// Authors:
//   Adam Strzelecki <adam.strzlecki@uj.edu.pl>
//   Piotr Bialas    <piotr.bialas@uj.edu.pl>
//
// Using Monte Carlo method and square detector scintilators.

#include <iostream>
#include <random>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "util/random.h"
#include "detector_ring.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"
#include "matrix_pixel_major.h"
#include "geometry/pixel.h"
#include "lor.h"
#include "model.h"
#include "util/png_writer.h"
#include "util/svg_ostream.h"
#include "util/util.h"

#include "monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/matrix.h"
#endif

// detect build variant
#if _OPENMP && HAVE_CUDA
#define VARIANT "OpenMP/CUDA"
#elif _OPENMP
#define VARIANT "OpenMP"
#elif HAVE_CUDA
#define VARIANT "CUDA"
#else
#define VARIANT "single-threaded CPU"
#endif

// all available detector shapes
typedef DetectorRing<double, int, SquareDetector<double>> SquareDetectorRing;
typedef DetectorRing<double, int, CircleDetector<double>> CircleDetectorRing;
typedef DetectorRing<double, int, TriangleDetector<double>>
    TriangleDetectorRing;
typedef DetectorRing<double, int, PolygonalDetector<6, double>>
    HexagonalDetectorRing;

template <typename DetectorRing, typename Model>
SparseMatrix<Pixel<>, LOR<>> run_cpu(cmdline::parser& cl,
                                     DetectorRing& detector_ring,
                                     Model& model);

template <typename DetectorRing>
void post_process(cmdline::parser& cl,
                  DetectorRing& detector_ring,
                  SparseMatrix<Pixel<>, LOR<>>& sparse_matrix);

void progress_callback(int pixel, int n_pixels);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    std::ostringstream msg;
    msg << "matrix_file ..." << std::endl;
    msg << "build: " << VARIANT << std::endl;
    msg << "note: All length options below should be expressed in meters.";
    cl.footer(msg.str());

    cl.add<cmdline::string>("config",
                            'c',
                            "load config file",
                            cmdline::dontsave,
                            cmdline::string(),
                            cmdline::default_reader<cmdline::string>(),
                            cmdline::load);
#if HAVE_CUDA
    cl.add("gpu", 'g', "run on GPU (via CUDA)");
    cl.add<int>("n-blocks", 0, "number of CUDA blocks", cmdline::dontsave, 64);
#endif
#if _OPENMP || HAVE_CUDA
    cl.add<int>(
        "n-threads", 't', "number of " VARIANT " threads", cmdline::dontsave);
#endif
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in one dimension", false, 256);
    cl.add<int>("m-pixel", 0, "starting pixel for partial matrix", false, 0);
    cl.add<int>("n-detectors", 'd', "number of detectors in ring", false);
    cl.add<int>("n-emissions",
                'e',
                "emissions per pixel",
                false,
                0,
                cmdline::default_reader<int>(),
                cmdline::not_from_file);
    cl.add<double>("radius", 'r', "inner detector ring radius", false);
    cl.add<double>("s-pixel", 'p', "pixel size", false);
    cl.add<double>(
        "tof-step", 'T', "TOF quantisation step for distance delta", false);
    cl.add<std::string>(
        "shape",
        'S',
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
        cmdline::oneof<std::string>(
            "always", "scintillator", /* obsolete */ "scintilator"));
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
    cl.add<tausworthe::seed_type>(
        "seed", 's', "random number generator seed", cmdline::dontsave);
    cl.add<cmdline::string>(
        "output",
        'o',
        "output binary triangular/full sparse system matrix",
        cmdline::dontsave);
    cl.add("full", 'f', "output full non-triangular sparse system matrix");

    // visual debugging params
    cl.add<cmdline::string>("png", 0, "output lor to png", cmdline::dontsave);
    cl.add<int>(
        "from", 0, "lor start detector to output", cmdline::dontsave, -1);
    cl.add<int>("to", 0, "lor end detector to output", cmdline::dontsave, -1);
    cl.add<int>("pos", 0, "position to output", cmdline::dontsave, -1);

    // printing & stats params
    cl.add("print", 0, "print triangular sparse system matrix");
    cl.add("stats", 0, "show stats");
    cl.add("wait", 0, "wait before exit");
    cl.add("verbose", 'v', "prints the iterations information on std::out");

    cl.try_parse(argc, argv);

    // convert obsolete acceptance option to length scale
    auto& length_scale = cl.get<double>("base-length");
    if (cl.exist("acceptance") && !cl.exist("base-length")) {
      length_scale = 1.0 / cl.get<double>("acceptance");
    }
    // FIXME: fixup for spelling mistake, present in previous versions
    auto& model = cl.get<std::string>("model");
    if (model == "scintilator") {
      model = "scintillator";
    }

    auto& n_pixels = cl.get<int>("n-pixels");
    auto& n_detectors = cl.get<int>("n-detectors");
    auto& n_emissions = cl.get<int>("n-emissions");
    auto& radius = cl.get<double>("radius");
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

    // load config files accompanying matrix files
    for (auto& fn : cl.rest()) {
      auto fn_sep = fn.find_last_of("\\/");
      auto fn_ext = fn.find_last_of(".");
      auto fn_wo_ext =
          fn.substr(0,
                    fn_ext != std::string::npos &&
                            (fn_sep == std::string::npos || fn_sep < fn_ext)
                        ? fn_ext
                        : std::string::npos);
      std::ifstream in(fn_wo_ext + ".cfg");
      if (!in.is_open())
        continue;
      // load except n-emissions
      auto n_prev_emissions = n_emissions;
      in >> cl;
      n_emissions = n_prev_emissions;
      break;  // only one config file allowed!
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

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
                 M_PI / std::atan2(d_detector, 2. * radius + d_detector / 2.)) /
             4) *
            4;
      } else {
        n_detectors =
            ((int)std::floor(M_PI / std::atan2(w_detector, 2. * radius)) / 4) *
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
  cl.exist("gpu") ? run_gpu(cl) : run_cpu(cl, detector_ring, model)
#else
#define _RUN(cl, detector_ring, model) run_cpu(cl, detector_ring, model)
#endif
#define RUN(detector_type, model_type, ...)                     \
  detector_type detector_ring(                                  \
      n_detectors, radius, w_detector, h_detector, d_detector); \
  model_type model{ __VA_ARGS__ };                              \
  auto sparse_matrix = _RUN(cl, detector_ring, model);          \
  post_process(cl, detector_ring, sparse_matrix)

    // run simmulation on given detector model & shape
    if (model == "always") {
      if (shape == "square") {
        RUN(SquareDetectorRing, AlwaysAccept<>);
      } else if (shape == "circle") {
        RUN(CircleDetectorRing, AlwaysAccept<>);
      } else if (shape == "triangle") {
        RUN(TriangleDetectorRing, AlwaysAccept<>);
      } else if (shape == "hexagon") {
        RUN(HexagonalDetectorRing, AlwaysAccept<>);
      }
    } else if (model == "scintillator") {
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
  }
  catch (cmdline::exception& ex) {
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
  }
  catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

template <typename DetectorRing, typename Model>
SparseMatrix<Pixel<>, LOR<>> run_cpu(cmdline::parser& cl,
                                     DetectorRing& detector_ring,
                                     Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& tof_step = cl.get<double>("tof-step");
  auto verbose = cl.exist("verbose");

  std::random_device rd;
  tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<tausworthe::seed_type>("seed"));
  }

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = detector_ring.n_positions(tof_step, max_bias);
  }

  typedef MatrixPixelMajor<Pixel<>, LOR<>> ComputeMatrix;
  ComputeMatrix::SparseMatrix sparse_matrix(
      n_pixels, n_detectors, n_tof_positions);

  for (auto& fn : cl.rest()) {
    ibstream in(fn, std::ios::binary);
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
        if (!cl.exist("n-detectors"))
          n_detectors = sparse_matrix.n_detectors();
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
    }

    catch (std::string& ex) {
      throw(ex + ": " + fn);
    }
    catch (const char* ex) {
      throw(std::string(ex) + ": " + fn);
    }
  }

  ComputeMatrix matrix(n_pixels, n_detectors, n_tof_positions);
  if (!sparse_matrix.empty()) {
    matrix << sparse_matrix;
    sparse_matrix.resize(0);
  }

  if (verbose) {
    std::cerr << "Monte-Carlo:" << std::endl;
#if _OPENMP
    std::cerr << " threads       = " << omp_get_max_threads() << std::endl;
#endif
    std::cerr << " pixels in row = " << n_pixels << std::endl;
    std::cerr << " outer radius  = " << detector_ring.outer_radius()
              << std::endl;
    std::cerr << " max bias      = " << max_bias << std::endl;
    std::cerr << " TOF step      = " << tof_step << std::endl;
    std::cerr << " TOF positions = " << n_tof_positions << std::endl;
    std::cerr << " emissions     = " << n_emissions << std::endl;
  }

#ifdef __linux__
  struct timespec start, stop;
  clock_gettime(CLOCK_REALTIME, &start);
#endif

  MonteCarlo<DetectorRing, ComputeMatrix> monte_carlo(
      detector_ring, matrix, s_pixel, tof_step, m_pixel);
  monte_carlo(gen, model, n_emissions, verbose ? progress_callback : NULL);

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
    auto fn = cl.get<cmdline::string>("output");
    auto fn_sep = fn.find_last_of("\\/");
    auto fn_ext = fn.find_last_of(".");
    auto fn_wo_ext =
        fn.substr(0,
                  fn_ext != std::string::npos &&
                          (fn_sep == std::string::npos || fn_sep < fn_ext)
                      ? fn_ext
                      : std::string::npos);
    auto fn_wo_path =
        fn_wo_ext.substr(fn_sep != std::string::npos ? fn_sep + 1 : 0);

    bool full = cl.exist("full");
    obstream out(fn, std::ios::binary | std::ios::trunc);
    if (full) {
      auto full_matrix = sparse_matrix.to_full();
      out << full_matrix;
    } else {
      out << sparse_matrix;
    }

    std::ofstream os(fn_wo_ext + ".cfg", std::ios::trunc);
    os << cl;

    try {
      png_writer png(fn_wo_ext + ".png");
      sparse_matrix.output_bitmap(png);
    }
    catch (const char* ex) {
      // don't bail out just produce warning
      std::cerr << "warning: " << ex << std::endl;
    }

    svg_ostream<> svg(fn_wo_ext + ".svg",
                      detector_ring.outer_radius(),
                      detector_ring.outer_radius(),
                      1024.,
                      1024.);
    svg << detector_ring;

    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2.,
                   -(s_pixel * n_pixels) / 2.,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);
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

    auto fn = cl.get<cmdline::string>("png");
    auto fn_sep = fn.find_last_of("\\/");
    auto fn_ext = fn.find_last_of(".");
    auto fn_wo_ext =
        fn.substr(0,
                  fn_ext != std::string::npos &&
                          (fn_sep == std::string::npos || fn_sep < fn_ext)
                      ? fn_ext
                      : std::string::npos);
    auto fn_wo_path =
        fn_wo_ext.substr(fn_sep != std::string::npos ? fn_sep + 1 : 0);

    png_writer png(fn);
    auto position = cl.get<int>("pos");
    if (cl.exist("full")) {
      sparse_matrix.to_full().output_bitmap(png, lor, position);
    } else {
      sparse_matrix.output_bitmap(png, lor, position);
    }

    svg_ostream<> svg(fn_wo_ext + ".svg",
                      detector_ring.outer_radius(),
                      detector_ring.outer_radius(),
                      1024.,
                      1024.);
    svg << detector_ring;

    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2.,
                   -(s_pixel * n_pixels) / 2.,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);
  }
}

void progress_callback(int pixel, int n_pixels) {
  static time_t start_time = 0;
  if (!start_time)
    start_time = time(NULL);
  report_progress(start_time, pixel, n_pixels);
}

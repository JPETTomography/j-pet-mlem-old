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
#include "matrix_pixel_major.h"
#include "geometry/pixel.h"
#include "lor.h"
#include "model.h"
#include "util/png_writer.h"
#include "util/svg_ostream.h"

#include "monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

#if HAVE_CUDA
#include "cuda/matrix.h"
#endif

// detect build variant
#if _OPENMP /* both */&& HAVE_CUDA
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

template <class DetectorRing> void run(cmdline::parser& cl);

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
    cl.add<int>("n-detectors", 'd', "number of detectors in ring", false, 64);
    cl.add<int>("n-emissions",
                'e',
                "emissions per pixel",
                false,
                0,
                cmdline::default_reader<int>(),
                cmdline::not_from_file);
    cl.add<double>("radius", 'r', "inner detector ring radius", false);
    cl.add<double>("s-pixel", 'p', "pixel size", false);
    cl.add<double>("tof-step", 'T', "TOF quantisation step", false);
    cl.add<std::string>(
        "shape",
        'S',
        "detector (scintillator) shape (square, circle,triangle)",
        false,
        "square",
        cmdline::oneof<std::string>("square", "circle", "triangle"));
    cl.add<double>("w-detector", 'w', "detector width", false);
    cl.add<double>("h-detector", 'h', "detector height", false);
    cl.add<std::string>("model",
                        'm',
                        "acceptance model (always, scintillator)",
                        false,
                        "scintillator",
                        cmdline::oneof<std::string>("always", "scintillator"));
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

#if HAVE_CUDA
    if (cl.exist("gpu")) {
      run_gpu(cl);
    } else
#endif
    {
      auto& shape = cl.get<std::string>("shape");

      // run simmulation on given detector shape
      if (shape == "square") {
        run<SquareDetectorRing>(cl);
      } else if (shape == "circle") {
        run<CircleDetectorRing>(cl);
      }
    }
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
  return 0;
}

template <class DetectorRing> void run(cmdline::parser& cl) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& n_detectors = cl.get<int>("n-detectors");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto& radius = cl.get<double>("radius");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& w_detector = cl.get<double>("w-detector");
  auto& h_detector = cl.get<double>("h-detector");
  auto& tof_step = cl.get<double>("tof-step");
  auto& model = cl.get<std::string>("model");
  auto& length_scale = cl.get<double>("base-length");
  auto verbose = cl.exist("verbose");

  // convert obsolete acceptance option to length scale
  if (cl.exist("acceptance") && !cl.exist("base-length")) {
    length_scale = 1.0 / cl.get<double>("acceptance");
  }
  // FIXME: fixup for spelling mistake, present in previous versions
  if (cl.exist("model") && model == "scintilator") {
    model = "scintillator";
  }
  // check options
  if (cl.exist("png") && !cl.exist("from")) {
    throw("need to specify output --png option when --from is specified");
  }
  if (!cl.exist("png") && cl.exist("from")) {
    throw("need to specify --from lor when output --png option is specified");
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

  // automatic pixel size
  if (!cl.exist("radius")) {
    if (!cl.exist("s-pixel")) {
      radius = M_SQRT2;  // exact result
    } else {
      radius = s_pixel * n_pixels / M_SQRT2;
    }
    std::cerr << "--radius=" << radius << std::endl;
  }

  // automatic radius
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radius")) {
      s_pixel = 2. / n_pixels;  // exact result
    } else {
      s_pixel = M_SQRT2 * radius / n_pixels;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radius / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::random_device rd;
  tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<tausworthe::seed_type>("seed"));
  }

  DetectorRing detector_ring(n_detectors, radius, w_detector, h_detector);

  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    if (model == "always")
      max_bias = AlwaysAccept<>::max_bias();
    else if (model == "scintilator")
      max_bias = ScintilatorAccept<>::max_bias();
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
      detector_ring, matrix, s_pixel, tof_step);
  if (model == "always")
    monte_carlo(gen, AlwaysAccept<>(), n_emissions);
  if (model == "scintillator")
    monte_carlo(gen, ScintilatorAccept<>(length_scale), n_emissions);

#ifdef __linux__
  if (verbose) {
    clock_gettime(CLOCK_REALTIME, &stop);
    std::cerr << "time : " << ((1.0e9 * stop.tv_sec + stop.tv_nsec) -
                               (1.0e9 * start.tv_sec + start.tv_nsec)) /
                                  1.0e9 << std::endl;
  }
#endif

  sparse_matrix = matrix.to_sparse();

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
    if (verbose) {
      std::cerr << "save sparse " << (full ? "full" : "triangular")
                << " matrix: " << fn << std::endl;
    }
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
      matrix.output_bitmap(png);
    }
    catch (const char* ex) {
      // don't bail out just produce warning
      std::cerr << "warning: " << ex << std::endl;
    }

    svg_ostream<> svg(fn_wo_ext + ".svg",
                      radius + h_detector,
                      radius + h_detector,
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
                      radius + h_detector,
                      radius + h_detector,
                      1024.,
                      1024.);
    svg << detector_ring;

    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2.,
                   -(s_pixel * n_pixels) / 2.,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);
  }

  // show stats if requested
  if (cl.exist("stats")) {
    auto pixel_max = 0;
    auto pixel_min = std::numeric_limits<decltype(pixel_max)>::max();
    for (auto y = 0; y < n_pixels; ++y) {
      for (auto x = 0; x < n_pixels; ++x) {
        auto hits = matrix[ComputeMatrix::Pixel(x, y)];
        pixel_min = std::min(pixel_min, hits);
        pixel_max = std::max(pixel_max, hits);
      }
    }
    std::cerr << "Non zero LORs: " << matrix.non_zero_lors() << '/'
              << detector_ring.lors() << std::endl;
    std::cerr << "Min hits: " << pixel_min << std::endl;
    std::cerr << "Max hits: " << pixel_max << std::endl;
  }

  if (cl.exist("print")) {
    sparse_matrix.sort_by_lor_n_pixel();
    std::cout << sparse_matrix;
  }

  if (cl.exist("wait")) {
    std::cerr << "Press Enter." << std::endl;
    while (getc(stdin) != '\n') {
    }
  }
}

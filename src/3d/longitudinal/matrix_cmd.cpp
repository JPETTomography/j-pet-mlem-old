/// \page cmd_3d_longitudinal_matrix 3d_longitudinal_matrix
/// \brief 3D Longitudinal PET system matrix construction tool
///
/// NYI
/// ===

#include <random>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "detector_set.h"
#include "2d/barrel/detector_set_builder.h"
#include "2d/barrel/detector_set.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/sparse_matrix.h"
#include "2d/barrel/matrix_pixel_major.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/model.h"
#include "monte_carlo.h"
#include "util/progress.h"

#include "util/png_writer.h"
#include "util/svg_ostream.h"

#include "options.h"
#if _OPENMP
#include <omp.h>
#endif

using namespace PET3D::Longitudinal;

using SquareScintillator = PET2D::Barrel::SquareDetector<float>;

using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareScintillator, 24, short>;
using LongitudinalDetectorSet = DetectorSet<DetectorSet2D>;

using Pixel = PET2D::Pixel<short>;
using LOR = PET2D::Barrel::LOR<short>;
using SparseMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR>;
using ComputeMatrix = PET2D::Barrel::MatrixPixelMajor<Pixel, LOR>;

template <typename DetectorRing, typename Model>
void print_parameters(cmdline::parser& cl, const DetectorRing& detector_ring);

template <typename Detector, typename Model>
static SparseMatrix run(cmdline::parser& cl,
                        Detector& detector_ring,
                        Model& model);

template <typename DetectorRing>
void post_process(cmdline::parser& cl,
                  DetectorRing& detector_ring,
                  SparseMatrix& sparse_matrix);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    add_matrix_options(cl);
    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

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
    calculate_detector_options(cl);

    const auto& model_name = cl.get<std::string>("model");
    const auto& length_scale = cl.get<double>("base-length");

    DetectorSet2D barrel =
        PET2D::Barrel::DetectorSetBuilder<DetectorSet2D>::buildMultipleRings(
            PET3D_LONGITUDINAL_DETECTOR_CL(cl, DetectorSet2D::F));
    LongitudinalDetectorSet detector_set = LongitudinalDetectorSet(barrel, 0.5);

    if (model_name == "always") {
      PET2D::Barrel::AlwaysAccept<float> model;
      auto sparse_matrix = run(cl, detector_set, model);
      post_process(cl, detector_set, sparse_matrix);
    } else if (model_name == "scintillator") {
      PET2D::Barrel::ScintillatorAccept<float> model(length_scale);
      auto sparse_matrix = run(cl, detector_set, model);
      post_process(cl, detector_set, sparse_matrix);
    } else {
      std::cerr << "unknown model : `" << model_name << "'\n";
      exit(-1);
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

template <typename Detector, typename Model>
static SparseMatrix run(cmdline::parser& cl,
                        Detector& detector_ring,
                        Model& model) {

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& m_pixel = cl.get<int>("m-pixel");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& n_emissions = cl.get<int>("n-emissions");
  auto verbose = cl.exist("verbose");

  std::random_device rd;
  util::random::tausworthe gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<util::random::tausworthe::seed_type>("seed"));
  }

  int n_tof_positions = 1;

  ComputeMatrix::SparseMatrix sparse_matrix(
      n_pixels, detector_ring.barrel.size(), n_tof_positions);

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

  ComputeMatrix matrix(n_pixels, detector_ring.barrel.size(), n_tof_positions);

#ifdef __linux__
  struct timespec start, stop;
  clock_gettime(CLOCK_REALTIME, &start);
#endif

  PET3D::Barrel::MonteCarlo<LongitudinalDetectorSet,
                            ComputeMatrix,
                            typename LongitudinalDetectorSet::F,
                            typename LongitudinalDetectorSet::S>
      monte_carlo(detector_ring, matrix, s_pixel, m_pixel);
  util::progress progress(verbose, matrix.total_n_pixels_in_triangle(), 1);
  monte_carlo(0.0f, gen, model, n_emissions, progress);

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

    util::svg_ostream<typename DetectorRing::F> svg(
        fn_wo_ext + ".svg",
        detector_ring.barrel.outer_radius(),
        detector_ring.barrel.outer_radius(),
        1024.,
        1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << const_cast<DetectorSet2D&>(detector_ring.barrel);
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
      sparse_matrix.to_full().output_bitmap(png, lor, position);
    } else {
      sparse_matrix.output_bitmap(png, lor, position);
    }

    util::svg_ostream<typename DetectorRing::F> svg(
        fn_wo_ext + ".svg",
        detector_ring.barrel.outer_radius(),
        detector_ring.barrel.outer_radius(),
        1024.,
        1024.);
    svg.link_image(fn_wo_path + ".png",
                   -(s_pixel * n_pixels) / 2,
                   -(s_pixel * n_pixels) / 2,
                   s_pixel * n_pixels,
                   s_pixel * n_pixels);

    svg << const_cast<DetectorSet2D&>(detector_ring.barrel);
  }
}

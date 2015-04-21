/// \page cmd_3d_longitudinal_matrix 3d_longitudinal_matrix
/// \brief 3D Longitudinal PET system matrix construction tool
///
/// NYI
/// ===

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

#include "options.h"
#if _OPENMP
#include <omp.h>
#endif

using namespace PET3D::Longitudinal;

using SquareScintillator = PET2D::Barrel::SquareDetector<float>;

using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareScintillator, 24>;
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
      run(cl, detector_set, model);
    } else if (model_name == "stintillator") {
      PET2D::Barrel::ScintillatorAccept<float> model(length_scale);
      run(cl, detector_set, model);
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

  int n_tof_positions = 1;

  ComputeMatrix matrix(n_pixels, detector_ring.barrel.size(), n_tof_positions);

  return matrix.to_sparse();
}

template <typename DetectorRing>
void post_process(cmdline::parser& cl,
                  DetectorRing& detector_ring,
                  SparseMatrix& sparse_matrix) {}

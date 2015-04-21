/// \page cmd_3d_longitudinal_matrix 3d_longitudinal_matrix
/// \brief 3D Longitudinal PET system matrix construction tool
///
/// NYI
/// ===

#include "detector_set.h"
#include "2d/barrel/detector_set.h"
#include "2d/barrel/square_detector.h"

#include "options.h"
#if _OPENMP
#include <omp.h>
#endif

using namespace PET3D::Longitudinal;

using SquareScintillator = PET2D::Barrel::SquareDetector<float>;

using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareScintillator, 24>;
using LongitudinalDetectorSet = DetectorSet<DetectorSet2D>;


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

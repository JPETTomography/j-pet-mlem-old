#include "2d/barrel/detector_set.h"
#include "2d/strip/detector.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET hybrid barrel & strip
namespace Hybrid {

/// Detector made of several other detectors
template <typename DetectorType = PET2D::Barrel::SquareDetector<double>,
          std::size_t MaxDetectors = MAX_DETECTORS,
          typename SType = int>
class DetectorSet
    : public PET2D::Barrel::DetectorSet<DetectorType, MaxDetectors, SType>,
      public PET2D::Strip::Detector<typename DetectorType::F> {
  // FIXME: implement me!
};
}
}

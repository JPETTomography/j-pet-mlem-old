#include "2d/barrel/detector_set.h"
#include "2d/strip/detector.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET hybrid barrel & strip
namespace Hybrid {

/// 3D detector made of several scintillators

/// This represents fusion of PET2D::Barrel:DetectorSet in `x-y` axis and
/// PET2D::Strip::Detector in `y-z` asis.
///
/// \image html detector3Daxis.pdf.png
///
/// \sa PET2D::Barrel:DetectorSet, PET2D::Strip::Detector

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

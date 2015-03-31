#include "2d/barrel/detector_set.h"

/// Three-dimensional PET
namespace PET3D {
/// Three-dimensional PET with longitudinal direction(z) added to the barrel
/// detector
namespace Longitudinal {

/// 3D detector made of several scintillators

template <typename DetectorType2D = PET2D::Barrel::SquareDetector<double>,
          std::size_t MaxDetectors = MAX_DETECTORS,
          typename FType = float,
          typename SType = int>
class DetectorSet {
  // FIXME: implement me!

  DetectorSet(const DetectorType2D& barrel_a, FType length_a)
      : barrel(barrel_a), length(length_a) {}

 private:
  const DetectorType2D& barrel;
  const FType length;
};

}  // Longitudinal
}  // PET3D

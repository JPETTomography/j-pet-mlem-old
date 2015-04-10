#ifndef DETECTORSETBUILDER_H
#define DETECTORSETBUILDER_H

#include "2d/barrel/detector_set.h"

namespace PET2D {
namespace Barrel {

template <typename DetectorType, std::size_t MaxDetectors, typename SType>
class DetectorSetBuilder {
 public:
  using FType = typename DetectorType::F;
  using DetectorSet =
      PET2D::Barrel::DetectorSet<DetectorType, MaxDetectors, SType>;

  static DetectorSet buildSingleRing(FType radius,
                                     SType n_detectors,
                                     FType w_detector,
                                     FType h_detector,
                                     FType d_detector = 0) {
    return DetectorSet(radius, n_detectors, w_detector, h_detector, d_detector);
  }
};
}
}

#endif  // DETECTORSETBUILDER_H

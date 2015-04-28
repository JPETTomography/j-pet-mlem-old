#pragma once

#include "2d/geometry/pixel.h"
namespace PET2D {
namespace Barrel {

enum class Axis { X, Y, XY };

template <typename SType> class SymmetryDescriptor {
 public:
  using S = SType;
  SymmetryDescriptor(int n_detectors, int n_symetries) {
    detectors_ = new S[n_detectors * n_symetries];
  }

  /**
   * @brief symmetric_detector
   * Returns symmetric detector on a ring of n_detectors, assuming that detector
   * zero is on the
   * positive X-axis (rotation=0).
   */
  S symmetric_detector(S detector, S symmetry) const {
    return detectors_[detector * symmetry + symmetry];
  };
  Pixel<S> pixel(Pixel<S> pixel, S symmetry);

  S ring_symmetric_detector(S n_detectors, S detector, S symmetry) const {
    if (symmetry & Axis::X) {
      detector = (n_detectors - detector) % n_detectors;  // x-axis
    }
    if (symmetry & Axis::Y) {
      detector =
          ((n_detectors + n_detectors / 2) - detector) % n_detectors;  // y-axis
    }
    if (symmetry & Axis::XY) {
      detector = ((n_detectors + n_detectors / 4) - detector) %
                 n_detectors;  // xy-axis
    }
    return detector;
  }

  /**
   * @brief symmetric_detector
   * Returns symmetric detector on a ring of n_detectors, assuming that detector
   * zero is
   * at the angle Pi/n_detectors (rotation = 0.5).
   */
  S rotated_ring_symmetric_detector(S n_detectors,
                                    S detector,
                                    S symmetry) const {
    if (symmetry & Axis::X) {
      detector = (n_detectors - (detector + 1)) % n_detectors;  // x-axis
    }
    if (symmetry & Axis::Y) {
      detector = ((n_detectors + n_detectors / 2) - (detector + 1)) %
                 n_detectors;  // y-axis
    }
    if (symmetry & Axis::XY) {
      detector = ((n_detectors + n_detectors / 4) - (detector + 1)) %
                 n_detectors;  // xy-axis
    }
    return detector;
  }

 private:
  S* detectors_;
};
}
}

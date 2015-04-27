#pragma once

#include "2d/geometry/pixel.h"
namespace PET2D {
namespace Barrel {

template <typename SType> class SymmetryDescriptor {
 public:
  using S = SType;
  SymmetryDescriptor(int n_detectors, int n_symetries) {
    detectors_ = new S[n_detectors * n_symetries]
  }

  S detector(S detector, S symmetry) const {
      return detectors_[detector*symmetry+ symmetry];

  };
  Pixel<S> pixel(Pixel<S> pixel, symmetry);

 private:
  S* detectors_;
}
}
}

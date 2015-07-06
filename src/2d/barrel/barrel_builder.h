#pragma once

#include "square_detector.h"
#include "scanner_builder.h"

namespace PET2D {
namespace Barrel {

using SquareDetectorType = PET2D::Barrel::SquareDetector<float>;

template <typename SType> struct BarrelBuilder {
  using Detector = SquareDetectorType;
  using S = SType;
  using SmallBarrel = PET2D::Barrel::GenericScanner<Detector, 192, S>;
  using BigBarrel = PET2D::Barrel::GenericScanner<Detector, 192, S>;

  static SmallBarrel make_small_barrel() {
    float width = 0.005;
    float height = 0.019;
    float r = 0.180 - height / 2;

    SmallBarrel barrel =
        PET2D::Barrel::ScannerBuilder<SmallBarrel>::build_single_ring(
            r, 24, width, height);
    barrel.set_fov_radius(0.150);

    return barrel;
  }

  static BigBarrel make_big_barrel() {
    float width = 0.007;
    float height = 0.019;
    float r1 = 0.430 - height / 2;
    float r2 = 0.475 - height / 2;
    float r3 = 0.575 - height / 2;

    BigBarrel barrel =
        PET2D::Barrel::ScannerBuilder<BigBarrel>::build_multiple_rings(
            { r1, r2, r3 }, { 0, 0.5f, 0.5 }, { 48, 48, 96 }, width, height);

    barrel.set_fov_radius(0.400);
    return barrel;
  }
};
}
}

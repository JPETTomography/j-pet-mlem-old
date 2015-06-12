#pragma once

#include "square_detector.h"
#include "scanner_builder.h"

namespace PET2D {
namespace Barrel {

using SquareDetectorType = PET2D::Barrel::SquareDetector<float>;
using SmallBarrelType =
    PET2D::Barrel::GenericScanner<SquareDetectorType, 192, int>;

using BigBarrelType =
    PET2D::Barrel::GenericScanner<SquareDetectorType, 192, int>;

inline SmallBarrelType buildSmallBarrel() {
  float width = 0.005;
  float height = 0.019;
  float r = 0.180 - height / 2;

  SmallBarrelType barrel =
      PET2D::Barrel::ScannerBuilder<SmallBarrelType>::build_single_ring(
          r, 24, width, height);
  barrel.set_fov_radius(0.150);

  return barrel;
}

inline BigBarrelType buildBigBarrel() {
  float width = 0.007;
  float height = 0.019;
  float r1 = 0.430 - height / 2;
  float r2 = 0.475 - height / 2;
  float r3 = 0.575 - height / 2;

  BigBarrelType barrel =
      PET2D::Barrel::ScannerBuilder<BigBarrelType>::build_multiple_rings(
          { r1, r2, r3 }, { 0, 0.5f, 0.5 }, { 48, 48, 96 }, width, height);

  barrel.set_fov_radius(0.400);
  return barrel;
}
}
}

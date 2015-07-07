#pragma once

#include "scanner_builder.h"

namespace PET2D {
namespace Barrel {

template <class DetectorClass, typename SType> struct BarrelBuilder {
  using Detector = DetectorClass;
  using F = typename Detector::F;
  using S = SType;
  using SmallBarrel = PET2D::Barrel::GenericScanner<Detector, 192, S>;
  using BigBarrel = PET2D::Barrel::GenericScanner<Detector, 192, S>;

  static SmallBarrel make_small_barrel() {
    F width = F(0.005);
    F height = F(0.019);
    F r = F(0.180) - height / 2;

    SmallBarrel barrel =
        PET2D::Barrel::ScannerBuilder<SmallBarrel>::build_single_ring(
            r, 24, width, height);
    barrel.set_fov_radius(F(0.150));

    return barrel;
  }

  static BigBarrel make_big_barrel() {
    F width = F(0.007);
    F height = F(0.019);
    F r1 = F(0.430) - height / 2;
    F r2 = F(0.475) - height / 2;
    F r3 = F(0.575) - height / 2;

    BigBarrel barrel =
        PET2D::Barrel::ScannerBuilder<BigBarrel>::build_multiple_rings(
            { r1, r2, r3 },
            { 0, F(0.5), F(0.5) },
            { 48, 48, 96 },
            width,
            height);

    barrel.set_fov_radius(F(0.400));
    return barrel;
  }
};
}
}

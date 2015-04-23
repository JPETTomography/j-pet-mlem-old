#pragma once

#include "2d/barrel/barrel_builder.h"
#include "detector_set.h"

namespace PET3D {
namespace Longitudinal {

  using SmallBarrelType2D = PET2D::Barrel::SmallBarrelType;

  using SmallBarrelType = DetectorSet<SmallBarrelType2D>;

  using BigBarrelType2D = PET2D::Barrel::BigBarrelType;
  using BigBarrelType =  DetectorSet<BigBarrelType2D>;

  SmallBarrelType buildSmallBarrel() {
    SmallBarrelType2D barrel = PET2D::Barrel::buildSmallBarrel();
    return SmallBarrelType(barrel, 0.300f);
  }

  BigBarrelType buildBigBarrel() {
    BigBarrelType2D barrel = PET2D::Barrel::buildBigBarrel();
    return BigBarrelType(barrel, 0.500f);
  }
}
}



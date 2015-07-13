#pragma once

namespace PET3D {

/// 3D voxel
template <typename SType> struct Voxel {
  using S = SType;
  Voxel(S ix, S iy, S iz) : ix(ix), iy(iy), iz(iz) {}

  const S ix;
  const S iy;
  const S iz;
};

}  // PET3D

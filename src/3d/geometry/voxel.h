#ifndef VOXEL
#define VOXEL

namespace PET3D {

template <typename SType> struct Voxel {
  using S = SType;
  Voxel(S ix, S iy, S iz) : ix(ix), iy(iy), iz(iz) {}

  const S ix;
  const S iy;
  const S iz;
};
}
#endif  // VOXEL

#ifndef SENSITIVITY_MAPPER
#define SENSITIVITY_MAPPER

#include "3d/geometry/voxel_set.h"

namespace PET3D {
namespace Hybrid {
template <typename Scanner> class SensitivityMapper {
 public:
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using VoxelSet = PET3D::VoxelSet<F, S>;

  SensitivityMapper(Scanner& scanner, VoxelSet& voxel_set, size_t n_emissions)
      : scanner(scanner), voxel_set(voxel_set), n_emissions(n_emissions){};

  template <typename RNG> void map(const Voxel<S>& voxel, RNG& rng) {


  }

  void map() {}

 private:
  Scanner& scanner;
  VoxelSet& voxel_set;
  size_t n_emissions;
};
}
}

#endif  // SENSITIVITY_MAPPER

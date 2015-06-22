#pragma once

#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_grid.h"

namespace PET3D {

template <typename FType, typename SType> class VoxelSet {
  using F = FType;
  using S = SType;
  using Voxel = PET3D::Voxel<S>;

 public:
  VoxelSet(const PET3D::VoxelGrid<F, S>& grid) : v_grid(grid) {}

  typename std::vector<Voxel>::const_iterator begin() const {
    return voxels_.begin();
  }
  typename std::vector<Voxel>::const_iterator end() const {
    return voxels_.end();
  }

  Voxel voxel(int i) const { return voxels_[i]; }
  F& value(int i) { return values_[i]; }
  size_t size() const { return voxels_.size(); }

  void push_back(const Voxel& voxel) {
    voxels_.push_back(voxel);
    values_.push_back(0);
  }

  const PET3D::VoxelGrid<F, S> v_grid;

 private:
  std::vector<Voxel> voxels_;
  std::vector<F> values_;
};
}

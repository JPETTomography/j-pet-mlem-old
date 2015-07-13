#pragma once

#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_grid.h"

namespace PET3D {

/// Set of 3D voxels
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

  void add_triangular_z_slice(S iz, F fov_radius) {
    auto& p_grid = v_grid.pixel_grid;
    auto cix = p_grid.n_columns / 2;
    auto ciy = p_grid.n_rows / 2;
    for (S ix = cix; ix < p_grid.n_columns; ix++)
      for (S iy = ciy; iy <= ix; iy++) {
        auto p = p_grid.center_at(ix, iy);
        if (p.x * p.x + p.y * p.y <= fov_radius * fov_radius) {
          push_back(Voxel(ix, iy, iz));
        }
      }
  }

  void add_y_slice(S iy, F fov_radius) {
    auto& p_grid = v_grid.pixel_grid;
    auto cix = p_grid.n_columns / 2;
    auto ciz = v_grid.n_planes / 2;
    for (S ix = cix; ix < p_grid.n_columns; ix++) {
      for (S iz = ciz; iz <= v_grid.n_planes; iz++) {
        auto p = p_grid.center_at(ix, iy);
        if (p.x * p.x + p.y * p.y <= fov_radius * fov_radius) {
          this->push_back(Voxel(ix, iy, iz));
        }
      }
    }
  }

  const PET3D::VoxelGrid<F, S> v_grid;

 private:
  std::vector<Voxel> voxels_;
  std::vector<F> values_;
};

}  // PET3D

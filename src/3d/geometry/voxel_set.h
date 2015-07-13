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
  VoxelSet(const PET3D::VoxelGrid<F, S>& grid) : grid(grid) {}

  typename std::vector<Voxel>::const_iterator begin() const {
    return voxels_.begin();
  }
  typename std::vector<Voxel>::const_iterator end() const {
    return voxels_.end();
  }

  Voxel voxel(int i) const { return voxels_[i]; }
  F& value(int i) { return values_[i]; }
  size_t size() const { return voxels_.size(); }

  template <typename... Args> void emplace_back(Args&&... args) {
    voxels_.emplace_back(std::forward<Args>(args)...);
    values_.push_back(0);
  }

  void push_back(const Voxel& voxel) {
    voxels_.push_back(voxel);
    values_.push_back(0);
  }

  void add_triangular_z_slice(S iz, F fov_radius) {
    auto& pixel_grid = grid.pixel_grid;
    auto cix = pixel_grid.n_columns / 2;
    auto ciy = pixel_grid.n_rows / 2;
    for (S ix = cix; ix < pixel_grid.n_columns; ix++)
      for (S iy = ciy; iy <= ix; iy++) {
        auto p = pixel_grid.center_at(ix, iy);
        if (p.x * p.x + p.y * p.y <= fov_radius * fov_radius) {
          emplace_back(ix, iy, iz);
        }
      }
  }

  void add_y_slice(S iy, F fov_radius) {
    auto& pixel_grid = grid.pixel_grid;
    auto cix = pixel_grid.n_columns / 2;
    auto ciz = grid.n_planes / 2;
    for (S ix = cix; ix < pixel_grid.n_columns; ix++) {
      for (S iz = ciz; iz <= grid.n_planes; iz++) {
        auto p = pixel_grid.center_at(ix, iy);
        if (p.x * p.x + p.y * p.y <= fov_radius * fov_radius) {
          emplace_back(ix, iy, iz);
        }
      }
    }
  }

  const PET3D::VoxelGrid<F, S> grid;

 private:
  std::vector<Voxel> voxels_;
  std::vector<F> values_;
};

}  // PET3D

#pragma once

#include "2d/geometry/pixel_grid.h"
#include "3d/geometry/point.h"

namespace PET3D {
template <typename FType, typename SType> class VoxelGrid {
  using F = FType;
  using S = SType;
  using PixelGrid = PET2D::PixelGrid<F, S>;
  using Point = PET3D::Point<F>;

 public:
  VoxelGrid(const PixelGrid& pixel_grid, F z_left, S n_planes)
      : pixel_grid(pixel_grid),
        z_left(z_left),
        n_planes(n_planes),
        n_voxels(pixel_grid.n_pixels * n_planes){};

  Point lower_left_at(S column, S row, S plane) const {

  };

  Point center_at(S column, S row, S plane) const {
    auto p2d = pixel_grid.center_at(column, row);
    F z = center_z_at(column, row, plane);
    return Point(p2d.x, p2d.y, z);
  };

  F center_z_at(S column, S row, S plane) const {
    return  (plane + F(0.5)) * pixel_grid.pixel_size + z_left;
  }


  int index (S column, S row, S plane) const {
    return pixel_grid.index(column, row)+plane*pixel_grid.n_pixels;
  }
  const PixelGrid pixel_grid;
  const F z_left;
  const S n_planes;
  const int n_voxels;
};
}

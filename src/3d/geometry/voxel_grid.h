#pragma once

#include "2d/geometry/pixel_grid.h"
#include "3d/geometry/point.h"
#include "3d/geometry/voxel.h"

namespace PET3D {

/// 3D voxel grid description
////
/// 3D voxel grid description, without actual voxels storage
template <typename FType, typename SType> class VoxelGrid {
 public:
  using F = FType;
  using S = SType;
  using PixelGrid = PET2D::PixelGrid<F, S>;
  using Pixel = typename PixelGrid::Pixel;
  using Point = PET3D::Point<F>;
  using Vector = PET3D::Vector<F>;
  using Voxel = PET3D::Voxel<S>;

  _ VoxelGrid(const PixelGrid& pixel_grid, F z_left, S n_planes)
      : pixel_grid(pixel_grid),
        z_left(z_left),
        n_planes(n_planes),
        n_voxels(pixel_grid.n_pixels * n_planes) {}

  _ Point lower_left_at(Voxel voxel) const {
    auto point2d = pixel_grid.lower_left_at(Pixel(voxel.x, voxel.y));
    F z = voxel.z * pixel_grid.pixel_size + z_left;
    return Point(point2d.x, point2d.y, z);
  }

  _ Point center_at(Voxel voxel) const {
    auto p2d = pixel_grid.center_at(Pixel(voxel.x, voxel.y));
    F z = center_z_at(voxel);
    return Point(p2d.x, p2d.y, z);
  }

  _ F center_z_at(Voxel voxel) const {
    return (voxel.z + F(0.5)) * pixel_grid.pixel_size + z_left;
  }

  _ Voxel voxel_at(Point p) const {
    auto lower_left =
        Point(pixel_grid.lower_left.x, pixel_grid.lower_left.y, z_left);
    auto v = p - lower_left;
    S column = static_cast<S>(compat::floor(v.x / pixel_grid.pixel_size));
    S row = static_cast<S>(compat::floor(v.y / pixel_grid.pixel_size));
    S plane = static_cast<S>(compat::floor(v.z / pixel_grid.pixel_size));
    return Voxel(column, row, plane);
  }

  _ int index(Voxel voxel) const {
    return pixel_grid.index(Pixel(voxel.x, voxel.y)) +
           voxel.z * pixel_grid.n_pixels;
    // Experiments with better "cache aware" indexing no real speed up:
    // return pixel_grid.index(column, row) * n_planes + plane;
  }

  _ bool contains(Voxel pixel) const {
    return pixel.x >= 0 && pixel.x < pixel_grid.n_columns &&  //
           pixel.y >= 0 && pixel.y < pixel_grid.n_rows &&     //
           pixel.z >= 0 && pixel.z < n_planes;                //
  }

  const PixelGrid pixel_grid;
  const F z_left;
  const S n_planes;
  const int n_voxels;
};

template <typename FType, typename SType>
class VariableVoxelSizeVoxelGrid : public VoxelGrid<FType, SType> {
 public:
  using F = FType;
  using S = SType;
  using Base = VoxelGrid<F, S>;
  using PixelGrid = typename Base::PixelGrid;
  using Pixel = typename Base::Pixel;
  using Voxel = typename Base::Voxel;
  using Point = typename Base::Point;
  _ VariableVoxelSizeVoxelGrid(const PixelGrid& pixel_grid,
                               F z_left,
                               S n_planes,
                               F voxel_size)
      : Base(pixel_grid, z_left, n_planes), voxel_size(voxel_size) {}

  _ Point lower_left_at(Voxel voxel) const {
    auto point2d = this->pixel_grid.lower_left_at(Pixel(voxel.x, voxel.y));
    F z = voxel.z * voxel_size + this->z_left;
    return Point(point2d.x, point2d.y, z);
  }

  _ Point center_at(Voxel voxel) const {
    auto p2d = this->pixel_grid.center_at(Pixel(voxel.x, voxel.y));
    F z = this->center_z_at(voxel);
    return Point(p2d.x, p2d.y, z);
  }

  _ F center_z_at(Voxel voxel) const {
    return (voxel.z + F(0.5)) * voxel_size + this->z_left;
  }

  _ Voxel voxel_at(Point p) const {
    auto lower_left = Point(this->pixel_grid.lower_left.x,
                            this->pixel_grid.lower_left.y,
                            this->z_left);
    auto v = p - lower_left;
    S column = static_cast<S>(compat::floor(v.x / this->pixel_grid.pixel_size));
    S row = static_cast<S>(compat::floor(v.y / this->pixel_grid.pixel_size));
    S plane = static_cast<S>(compat::floor(v.z / voxel_size));
    return Voxel(column, row, plane);
  }

  const F voxel_size;
};

}  // PET3D

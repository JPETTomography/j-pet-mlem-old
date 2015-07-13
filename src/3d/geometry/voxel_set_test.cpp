#include "util/test.h"

#include "3d/geometry/voxel_set.h"

#include "common/mathematica_graphics.h"

using F = float;
using S = int;
using MathematicaGraphics = Common::MathematicaGraphics<F>;

TEST("3d/geometry/voxel_set") {
  using Voxel = PET3D::Voxel<S>;
  using VoxelSet = PET3D::VoxelSet<F, S>;
  using VoxelGrid = PET3D::VoxelGrid<F, S>;
  using Point = PET2D::Point<F>;
  using PixelGrid = PET2D::PixelGrid<F, S>;

  PixelGrid pixel_grid(80, 80, 0.005, Point(-0.200, -0.200));
  VoxelGrid voxel_grid(pixel_grid, -0.200, 80);

  SECTION("empty") {
    VoxelSet voxel_set(voxel_grid);
    REQUIRE(voxel_set.size() == 0);
  }

  SECTION("push_back") {
    VoxelSet voxel_set(voxel_grid);
    voxel_set.push_back(Voxel(1, 2, 3));
    REQUIRE(voxel_set.size() == 1);
    for (auto& voxel : voxel_set) {
      REQUIRE(voxel.ix == 1);
      REQUIRE(voxel.iy == 2);
      REQUIRE(voxel.iz == 3);
    }
  }

  SECTION("triangular_z_slice") {
    VoxelSet voxel_set(voxel_grid);
    voxel_set.add_triangular_z_slice(41, 0.200);
    std::ofstream out("test_output/triangular_voxels.m");
    MathematicaGraphics graphics(out);
    for (auto& voxel : voxel_set) {
      graphics.add_pixel(pixel_grid, voxel.ix, voxel.iy);
    }
  }

  SECTION("y_slice") {
    VoxelSet voxel_set(voxel_grid);
    voxel_set.add_y_slice(79, 0.200);
    std::ofstream out("test_output/yslice_voxels.m");
    MathematicaGraphics graphics(out);
    for (auto& voxel : voxel_set) {
      graphics.add_pixel(pixel_grid, voxel.iz, voxel.ix);
    }
  }
}

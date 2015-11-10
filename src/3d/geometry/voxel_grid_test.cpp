#include "util/test.h"

#include "3d/geometry/voxel_grid.h"

#include "common/types.h"

using Voxel = PET3D::Voxel<S>;

TEST("3d/geometry/voxel_grid") {
  PET2D::PixelGrid<F, S> pixel_grid(10, 8, 0.005, PET2D::Point<F>(0, 0));
  PET3D::VoxelGrid<F, S> grid(pixel_grid, -0.015, 6);
  REQUIRE(grid.n_voxels == 10 * 8 * 6);
  {
    auto p = grid.center_at(Voxel(1, 2, 3));
    REQUIRE(p.x == 0.0075_e7);
    REQUIRE(p.y == 0.0125_e7);
    REQUIRE(p.z == 0.0025_e7);
  }
  {
    auto p = grid.lower_left_at(Voxel(1, 2, 3));
    REQUIRE(p.x == 0.005_e7);
    REQUIRE(p.y == 0.010_e7);
    REQUIRE(p.z == 0.000_e7);
  }
}

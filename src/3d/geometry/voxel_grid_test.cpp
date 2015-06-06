#include "util/test.h"

#include "3d/geometry/voxel_grid.h"

TEST("3d/voxel_grid") {

  PET2D::PixelGrid<float, short> pixel_grid(
      10, 8, 0.005, PET2D::Point<float>(0, 0));
  PET3D::VoxelGrid<float, short> grid(pixel_grid, -0.015, 6);

  REQUIRE(grid.n_voxels == 10 * 8 * 6);
  {
    auto p = grid.center_at(1, 2, 3);
    REQUIRE(p.x == 0.0075_e7);
    REQUIRE(p.y == 0.0125_e7);
    REQUIRE(p.z == 0.0025_e7);
  }

  {
    auto p = grid.lower_left_at(1, 2, 3);
    REQUIRE(p.x == 0.005_e7);
    REQUIRE(p.y == 0.010_e7);
    REQUIRE(p.z == 0.000_e7);
  }
}

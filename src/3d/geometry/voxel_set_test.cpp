#include "util/test.h"

#include "3d/geometry/voxel_set.h"

using F = float;
using S = int;

TEST("VoxelSet") {
  using Voxel = PET3D::Voxel<S>;

  PET2D::PixelGrid<F, S> p_grid(80, 80, 0.005, PET2D::Point<F>(-0.200, 0.200));
  PET3D::VoxelGrid<F, S> v_grid(p_grid, -0.200, 80);

  VoxelSet<F, S> voxel_set(v_grid);

  REQUIRE(voxel_set.size() == 0);

  SECTION("push_back") {
    voxel_set.push_back(Voxel(1, 2, 3));
    REQUIRE(voxel_set.size() == 1);
    for (auto& voxel : voxel_set) {
      REQUIRE(voxel.ix == 1);
      REQUIRE(voxel.iy == 2);
      REQUIRE(voxel.iz == 3);
    }
  }
}

#include "util/test.h"

#include "3d/geometry/voxel_set.h"
#include "3d/geometry/voxel_set_builder.h"
#include "util/grapher.h"

using F = float;
using S = int;

TEST("3d/geometry/voxel_set") {
  using Voxel = PET3D::Voxel<S>;

  PET2D::PixelGrid<F, S> p_grid(80, 80, 0.005, PET2D::Point<F>(-0.200, -0.200));
  PET3D::VoxelGrid<F, S> v_grid(p_grid, -0.200, 80);

  PET3D::VoxelSet<F, S> voxel_set(v_grid);

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

  SECTION("BuildTriagularZSlice") {

    PET3D::VoxelSetBuilder<F, S>::BuildTriagularZSlice(voxel_set, 41, 0.200);

    std::ofstream out("test_output/triangular_voxels.m");
    Graphics<F> graphics(out);
    for (auto& voxel : voxel_set) {
      graphics.add_pixel(p_grid, voxel.ix, voxel.iy);
    }
  }

  SECTION("BuildYSlice") {
    PET3D::VoxelSetBuilder<F, S>::BuildYSlice(voxel_set, 79, 0.200);
    std::ofstream out("test_output/yslice_voxels.m");
    Graphics<F> graphics(out);
    for (auto& voxel : voxel_set) {
      graphics.add_pixel(p_grid, voxel.iz, voxel.ix);
    }
  }
}

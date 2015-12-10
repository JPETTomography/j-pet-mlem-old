#include <cmath>

#include "util/test.h"

#include "common/types.h"

#include "psf.h"

using Voxel = PET3D::Voxel<S>;
using VoxelMap = PET3D::VoxelMap<Voxel, F>;
using Vector = PET3D::Vector<F>;

static F data_3x3_point[] = {
  0, 0, 0,  //
  0, 0, 0,  //
  0, 0, 0,  //

  0, 0, 0,  //
  0, 2, 0,  //
  0, 0, 0,  //

  0, 0, 0,  //
  0, 0, 0,  //
  0, 0, 0,  //
};

static F data_3x3_fuzz[] = {
  0, 0, 0,  //
  0, 1, 0,  //
  0, 0, 0,  //

  0, 1, 0,  //
  1, 2, 1,  //
  0, 1, 0,  //

  0, 0, 0,  //
  0, 1, 0,  //
  0, 0, 0,  //
};

static F data_4x4_point_2x2[] = {
  1, 1, 1, 1,  //
  1, 1, 1, 1,  //
  1, 1, 1, 1,  //
  1, 1, 1, 1,  //

  1, 1, 1, 1,  //
  1, 6, 6, 1,  //
  1, 6, 6, 1,  //
  1, 1, 1, 1,  //

  1, 1, 1, 1,  //
  1, 6, 6, 1,  //
  1, 6, 6, 1,  //
  1, 1, 1, 1,  //

  1, 1, 1, 1,  //
  1, 1, 1, 1,  //
  1, 1, 1, 1,  //
  1, 1, 1, 1,  //
};

static F data_5x5_fuzz[] = {
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //

  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 1, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //

  0, 0, 0, 0, 0,  //
  0, 0, 1, 0, 0,  //
  0, 1, 2, 1, 0,  //
  0, 0, 1, 0, 0,  //
  0, 0, 0, 0, 0,  //

  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 1, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //

  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
  0, 0, 0, 0, 0,  //
};

TEST("3d/hybrid/psf/max") {
  F max;
  Voxel max_voxel;
  {
    VoxelMap map(3, 3, 3, data_3x3_point);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    CHECK(max == 2);
    CHECK(max_voxel == Voxel(1, 1, 1));
  }
  {
    VoxelMap map(3, 3, 3, data_3x3_fuzz);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    CHECK(max == 2);
    CHECK(max_voxel == Voxel(1, 1, 1));
  }
  {
    VoxelMap map(4, 4, 4, data_4x4_point_2x2);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    CHECK(max == 6);
    CHECK(max_voxel == Voxel(1, 1, 1));
  }
}

TEST("3d/hybrid/psf/left_right_above") {
  F max;
  Voxel max_voxel;
  Voxel left_above, right_above;
  {
    VoxelMap map(3, 3, 3, data_3x3_point);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    CHECK(left_above == Voxel(1, 1, 1));
    CHECK(right_above == Voxel(1, 1, 1));
  }
  {
    VoxelMap map(3, 3, 3, data_3x3_fuzz);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    CHECK(left_above == Voxel(0, 0, 0));
    CHECK(right_above == Voxel(2, 2, 2));
  }
  {
    VoxelMap map(4, 4, 4, data_4x4_point_2x2);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    CHECK(left_above == Voxel(1, 1, 1));
    CHECK(right_above == Voxel(2, 2, 2));
  }
}

TEST("3d/hybrid/psf/psf") {
  F max;
  Voxel max_voxel;
  Voxel left_above, right_above;
  Vector left, right, psf;
  {
    VoxelMap map(3, 3, 3, data_3x3_point);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    PET3D::Hybrid::PSF::calculate(
        map, max_voxel, max, left_above, right_above, left, right, psf);
    CHECK((left - Vector(.5, .5, .5)).length2() == 0._e7);
    CHECK((right - Vector(1.5, 1.5, 1.5)).length2() == 0._e7);
    CHECK((psf - Vector(1, 1, 1)).length2() == 0._e7);
  }
  {
    VoxelMap map(3, 3, 3, data_3x3_fuzz);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    PET3D::Hybrid::PSF::calculate(
        map, max_voxel, max, left_above, right_above, left, right, psf);
    CHECK((left - Vector(0, 0, 0)).length2() == 0._e7);
    CHECK((right - Vector(2, 2, 2)).length2() == 0._e7);
    CHECK(psf == Vector(-1, -1, -1));
  }
  {
    VoxelMap map(4, 4, 4, data_4x4_point_2x2);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    PET3D::Hybrid::PSF::calculate(
        map, max_voxel, max, left_above, right_above, left, right, psf);
    CHECK((left - Vector(0.4, 0.4, 0.4)).length2() == 0._e7);
    CHECK((right - Vector(2.6, 2.6, 2.6)).length2() == 0._e7);
    CHECK((psf - Vector(2.2, 2.2, 2.2)).length2() == 0._e7);
  }
  {
    VoxelMap map(5, 5, 5, data_5x5_fuzz);
    PET3D::Hybrid::PSF::find_max(map, max_voxel, max);
    PET3D::Hybrid::PSF::find_left_right_above_half(
        map, max_voxel, max, left_above, right_above);
    PET3D::Hybrid::PSF::calculate(
        map, max_voxel, max, left_above, right_above, left, right, psf);
    CHECK((left - Vector(1, 1, 1)).length2() == 0._e7);
    CHECK((right - Vector(3, 3, 3)).length2() == 0._e7);
    CHECK(psf == Vector(2, 2, 2));
  }
}

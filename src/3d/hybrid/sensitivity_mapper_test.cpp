#include "util/test.h"

#include "3d/hybrid/sensitivity_mapper.h"
#include "3d/hybrid/scanner.h"
#include "2d/barrel/barrel_builder.h"
#include "util/random.h"
#include "2d/barrel/model.h"

using F = float;
using S = int;

TEST("SensitivityMapper") {
  using Voxel = PET3D::Voxel<S>;
  using BarrelType = PET2D::Barrel::BigBarrelType;
  using Scanner = PET3D::Hybrid::Scanner<BarrelType>;

  PET2D::PixelGrid<F, S> p_grid(80, 80, 0.005, PET2D::Point<F>(-0.200, -0.200));
  PET3D::VoxelGrid<F, S> v_grid(p_grid, -0.200, 80);

  PET3D::VoxelSet<F, S> voxel_set(v_grid);

  REQUIRE(voxel_set.size() == 0);
  voxel_set.push_back(Voxel(41, 41, 41));

  auto barrel = PET2D::Barrel::buildBigBarrel();
  Scanner scanner(barrel, 0.500);

  PET3D::Hybrid::SensitivityMapper<Scanner> mapper(scanner, voxel_set);

  PET2D::Barrel::ScintillatorAccept<F> scintillator(0.100);

  util::random::tausworthe gen(12212);
  mapper.map(gen, scintillator, 10000);

  for (int i = 0; i < voxel_set.size(); ++i) {
    std::cout << voxel_set.value(i) << "\n";
  }
}

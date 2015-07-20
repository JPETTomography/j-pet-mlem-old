#include "util/test.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/options.h"
#include "3d/hybrid/sensitivity_mapper.h"
#include "3d/hybrid/scanner.h"
#include "util/random.h"

#include "common/model.h"
#include "common/types.h"

using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner2DBuilder = PET2D::Barrel::ScannerBuilder<Scanner2D>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;

static Scanner2D make_big_barrel() {
  std::ifstream in_cfg("samples/config/big.cfg");
  cmdline::parser cl;
  PET2D::Barrel::add_scanner_options(cl);
  cl.parse(in_cfg);
  return Scanner2DBuilder::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, Detector::F));
}

TEST("3d/hybrid/sensitivity_mapper") {
  using Voxel = PET3D::Voxel<S>;

  PET2D::PixelGrid<F, S> p_grid(80, 80, 0.005, PET2D::Point<F>(-0.200, -0.200));
  PET3D::VoxelGrid<F, S> v_grid(p_grid, -0.200, 80);

  PET3D::VoxelSet<F, S> voxel_set(v_grid);

  REQUIRE(voxel_set.size() == 0);
  voxel_set.push_back(Voxel(41, 41, 41));

  auto barrel = make_big_barrel();
  Scanner scanner(barrel, 0.500);

  PET3D::Hybrid::SensitivityMapper<Scanner> mapper(scanner, voxel_set);

  Common::ScintillatorAccept<F> scintillator(0.100);

  util::random::tausworthe gen(12212);
  mapper.map(gen, scintillator, 10000);
}

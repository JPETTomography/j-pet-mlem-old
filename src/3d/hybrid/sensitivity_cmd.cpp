#if _OPENMP
#include <omp.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/json.h"
#include "util/backtrace.h"

#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/options.h"
#include "3d/hybrid/scanner.h"
#include "3d/hybrid/sensitivity_mapper.h"
#include "util/random.h"
#include "common/model.h"

using F = float;
using S = int;

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Scanner2DBuilder = PET2D::Barrel::ScannerBuilder<Scanner2D>;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET2D::Barrel::add_scanner_options(cl);
  cl.add<int>("z-plane", 0, "z plane trianguler cut", false);
  cl.add<int>("y-plane", 0, "y plane cut", false);
  cl.add<int>(
      "n-pixels", 'n', "number of pixels in x and y  directions", false, 80);
  cl.add<int>("n-planes", '\0', "number pf z planes", false, 80);
  cl.add<float>("s-pixel", 'p', "voxel size", false, 0.005);
  cl.add<float>("fov-radius", '\0', "field of view radius", false, 0.400);
  cl.add<int>("n-emissions", 'e', "number of emission", false, 0);
  cl.add<cmdline::path>(
      "output", 'o', "output files template", false, "out.bin");
  cl.add<int>("n-threads", 'T', "number of threads", false);

  cl.try_parse(argc, argv);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto barrel = Scanner2DBuilder::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, SquareDetector::F));
  Scanner scanner(barrel, 0.500);
  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<float>("s-pixel");
  float ll = -s_pixel * n_pixels / 2;
  PET2D::PixelGrid<F, S> pixel_grid(
      n_pixels, n_pixels, s_pixel, PET2D::Point<F>(ll, ll));

  int n_planes = cl.get<int>("n-planes");
  PET3D::VoxelGrid<F, S> voxel_grid(
      pixel_grid, -s_pixel * n_planes / 2, n_planes);

  PET3D::VoxelSet<F, S> voxel_set(voxel_grid);
  if (cl.exist("z-plane")) {
    voxel_set.add_triangular_z_slice(cl.get<int>("z-plane"),
                                     cl.get<float>("fov-radius"));
  } else if (cl.exist("y-plane")) {
    voxel_set.add_y_slice(cl.get<int>("y-plane"), cl.get<float>("fov-radius"));
  } else {
    throw("you must specify --y-plane or --z-plane");
  }

  PET3D::Hybrid::SensitivityMapper<Scanner> mapper(scanner, voxel_set);

  Common::ScintillatorAccept<F> scintillator(0.100);

  util::random::tausworthe gen(12212);
  mapper.map(gen, scintillator, cl.get<int>("n-emissions"));

  if (cl.exist("output")) {
    auto fn = cl.get<cmdline::path>("output");
    auto fn_wo_ext = fn.wo_ext();

    util::obstream out(fn);
    for (size_t i = 0; i < voxel_set.size(); ++i) {
      out << voxel_set.voxel(i) << voxel_set.value(i);
    }

    std::ofstream out_json(fn_wo_ext + ".json");
    out_json << json(scanner.barrel);
  }

  CMDLINE_CATCH
}

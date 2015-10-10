/// \page cmd_3d_hybrid_sensitivity 3d_hybrid_sensitivity
/// \brief 3D Longitudinal PET sensitivity map construction tool
///
/// Creates sensitivity map for given scanner geometry.
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/matrix_cmd.txt
///
/// \sa \ref cmd_3d_hybrid_phantom, \ref cmd_3d_hybrid_matrix

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

#include "scanner.h"
#include "sensitivity_mapper.h"
#include "options.h"

#include "util/random.h"

#include "common/model.h"
#include "common/types.h"

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Scanner2DBuilder = PET2D::Barrel::ScannerBuilder<Scanner2D>;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  PET3D::Hybrid::add_sensitivity_options(cl);
  cl.parse_check(argc, argv);
  PET3D::Hybrid::calculate_scanner_options(cl, argc);

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  Scanner scanner(
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_multiple_rings(
          PET3D_LONGITUDINAL_SCANNER_CL(cl, F)),
      F(cl.get<double>("length")));

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");
  auto ll = -s_pixel * n_pixels / 2;
  PET2D::PixelGrid<F, S> pixel_grid(
      n_pixels, n_pixels, s_pixel, PET2D::Point<F>(ll, ll));

  int n_planes = cl.get<int>("n-planes");
  PET3D::VoxelGrid<F, S> voxel_grid(
      pixel_grid, -s_pixel * n_planes / 2, n_planes);

  PET3D::VoxelSet<F, S> voxel_set(voxel_grid);
  if (cl.exist("z-plane")) {
    voxel_set.add_triangular_z_slice(cl.get<int>("z-plane"),
                                     cl.get<double>("fov-radius"));
  } else if (cl.exist("y-plane")) {
    voxel_set.add_y_slice(cl.get<int>("y-plane"), cl.get<double>("fov-radius"));
  } else {
    throw("you must specify --y-plane or --z-plane");
  }

  PET3D::Hybrid::SensitivityMapper<Scanner> mapper(scanner, voxel_set);

  Common::ScintillatorAccept<F> scintillator(F(0.100));

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

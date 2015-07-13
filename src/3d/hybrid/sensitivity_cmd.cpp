#if _OPENMP
#include <omp.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/json_ostream.h"

#include "2d/barrel/barrel_builder.h"
#include "3d/hybrid/scanner.h"
#include "3d/hybrid/sensitivity_mapper.h"
#include "util/random.h"
#include "common/model.h"

using F = float;
using S = int;

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using BarrelBuilder = PET2D::Barrel::BarrelBuilder<SquareDetector, S>;
using Scanner2D = BarrelBuilder::BigBarrel;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    cl.add<int>("z-plane", 0, "z plane trianguler cut", false);
    cl.add<int>("y-plane", 0, "y plane cut", false);
    cl.add<int>(
        "n-pixels", 'n', "number of pixels in x and y  directions", false, 80);
    cl.add<int>("n-planes", '\0', "number pf z planes", false, 80);
    cl.add<float>("pixel-size", 'p', "voxel size", false, 0.005);
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

    Scanner2D barrel = BarrelBuilder::make_big_barrel();
    Scanner scanner(barrel, 0.500);
    auto n_pixels = cl.get<int>("n-pixels");
    auto pixel_size = cl.get<float>("pixel-size");
    float ll = -pixel_size * n_pixels / 2;
    PET2D::PixelGrid<F, S> p_grid(
        n_pixels, n_pixels, pixel_size, PET2D::Point<F>(ll, ll));
    int n_planes = cl.get<int>("n-planes");
    PET3D::VoxelGrid<F, S> v_grid(p_grid, -pixel_size * n_planes / 2, n_planes);

    PET3D::VoxelSet<F, S> voxel_set(v_grid);
    if (cl.exist("z-plane")) {
      voxel_set.add_triangular_z_slice(cl.get<int>("z-plane"),
                                       cl.get<float>("fov-radius"));
    } else if (cl.exist("y-plane")) {
      voxel_set.add_y_slice(cl.get<int>("y-plane"),
                            cl.get<float>("fov-radius"));
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

      util::json_ostream json(fn_wo_ext + ".json");
      json << scanner.barrel;
    }

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

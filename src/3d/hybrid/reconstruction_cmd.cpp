#include <iostream>
#include <fstream>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "2d/barrel/options.h"
#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/lors_pixels_info.h"
#include "2d/strip/gausian_kernel.h"
#include "3d/hybrid/scanner.h"
#include "3d/hybrid/reconstruction.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

using F = float;
using S = int;
using RNG = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 192, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Point = PET2D::Point<F>;

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    cl.add<std::string>("lor-info", 0, "lor-pixel information", true);
    cl.add<float>("sigma-z", 0, "sigma-z", false, 0.015);
    cl.add<float>("sigma-dl", 0, "sigma-dl", false, 0.060);

    cl.add<double>("length", 0, "length of the detector", false, 0.3);
    cl.add<std::string>("response", 0, "detector responses", true);

    PET2D::Barrel::add_matrix_options(cl);
    cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
    cl.add<int>(
        "iterations", 'I', "number of iterations (per block)", false, 1);
#if _OPENMP
    cl.add<int>("n-threads",
                'T',
                "number of OpenMP threads to use",
                false,
                omp_get_max_threads());
#endif
    cl.try_parse(argc, argv);

    PET3D::Hybrid::set_big_barrel_options(cl);

    Scanner scanner = Scanner::build_scanner_from_cl(cl);
    scanner.set_sigmas(cl.get<float>("sigma-z"), cl.get<float>("sigma-dl"));
    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    auto lor_info_file_name = cl.get<std::string>("lor-info");
    int n_detectors;
    std::ifstream lor_info_istream(lor_info_file_name, std::ios::binary);
    lor_info_istream.read((char*)&n_detectors, sizeof(n_detectors));

    auto grid = PET2D::PixelGrid<F, S>::read(lor_info_istream);

    std::cout << grid.n_columns << "x" << grid.n_rows << " " << grid.pixel_size
              << "\n";
    if (n_detectors != scanner.barrel.size()) {
      throw("n_detectors mismatch");
    }
    PET2D::Barrel::LORsPixelsInfo<F, S> lor_info(scanner.barrel.size(), grid);
    lor_info.read(lor_info_istream);

    Reconstruction<Scanner, PET2D::Strip::GaussianKernel<F>> reconstruction(
        scanner, lor_info, -0.200, 80);

    std::ifstream response_stream(cl.get<std::string>("response"));
    reconstruction.fscanf_responses(response_stream);

    {
      std::ofstream gout("event.m");
      Graphics<float> graphics(gout);
      graphics.add(scanner.barrel);
      reconstruction.graph_frame_event(graphics, 0);
    }

    auto n_blocks = cl.get<int>("blocks");
    auto n_iter = cl.get<int>("iterations");

    for (int block = 0; block < n_blocks; ++block) {
      for (int i = 0; i < n_iter; i++) {
        std::cout << block * n_iter + i << " " << reconstruction.iterate()
                  << "\n";
      }
      char rho_file_name[64];
      sprintf(rho_file_name,
              "%s_%03d.bin",
              output_base_name.c_str(),
              (block + 1) * n_iter);
      std::ofstream out(rho_file_name);
      out.write((char*)&(*reconstruction.rho_begin()),
                reconstruction.n_voxels * sizeof(F));
    }
    std::cout << reconstruction.event_count() << " "
              << reconstruction.voxel_count() << " "
              << reconstruction.pixel_count() << "\n";
    std::cout << (double)reconstruction.voxel_count() /
                     reconstruction.event_count()
              << " ";
    std::cout << (double)reconstruction.pixel_count() /
                     reconstruction.event_count()
              << "\n";

    std::ofstream out("rho.bin");
    out.write((char*)&(*reconstruction.rho_begin()),
              reconstruction.n_voxels * sizeof(F));

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

  return 0;
}

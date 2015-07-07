#include <iostream>
#include <fstream>
#include <random>

#include "cmdline.h"

#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/bstream.h"

#include "2d/barrel/barrel_builder.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/options.h"
#include "2d/barrel/lors_pixels_info.h"
#include "2d/barrel/sparse_matrix.h"
#include "2d/barrel/lm_reconstruction.h"

#include "util/grapher.h"

#if _OPENMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_thread_num() 0
#endif

using F = float;
using S = int16_t;
using H = int;
using RNG = std::mt19937;
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using BarrelBuilder = PET2D::Barrel::BarrelBuilder<SquareDetector, S>;

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    cl.add<cmdline::path>("lor-info", 0, "lor-pixel information", true);
    cl.add<cmdline::path>("system", 0, "system maxtrix", false);
    cl.add<double>("sigma", 0, "sigma dl", false, 0.060);

    cl.add<double>("length", 0, "length of the detector", false, 0.3);
    cl.add<std::string>("response", 0, "detector responses", true);

    PET2D::Barrel::add_matrix_options(cl);
    cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
    cl.add<int>(
        "iterations", 'I', "number of iterations (per block)", false, 1);
    cl.add("graph", 'g', "make a graph", false);
    cl.add("event", 0, "event number", false, 0);
#if _OPENMP

#endif
    cl.try_parse(argc, argv);

    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    bool verbose = cl.exist("verbose");

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    S n_detectors;

    auto lor_info_file_name = cl.get<cmdline::path>("lor-info");
    std::ifstream lor_info_istream(lor_info_file_name, std::ios::binary);
    lor_info_istream.read((char*)&n_detectors, sizeof(n_detectors));

    std::cout << n_detectors << "\n";

    auto grid = PET2D::PixelGrid<F, S>::read(lor_info_istream);

    if (verbose)
      std::cout << grid.n_columns << "x" << grid.n_rows << " "
                << grid.pixel_size << "\n";

    PET2D::Barrel::LORsPixelsInfo<F, S> lor_info(n_detectors, grid);
    lor_info.read(lor_info_istream);
    if (verbose)
      std::cout << "read in lor_info\n";

    if (!cl.exist("system")) {

    } else {
      lor_info.erase_pixel_info();
      auto system_matrix_file_name = cl.get<cmdline::path>("system");
      util::ibstream system_matrix_istream(system_matrix_file_name);
      PET2D::Barrel::SparseMatrix<Pixel, LOR, S, H> matrix(
          system_matrix_istream);
      if (verbose)
        std::cout << "read in system matrix" << std::endl;
      matrix.sort_by_lor_n_pixel();
      matrix.merge_duplicates();
      F n_emissions = F(matrix.n_emissions());
      if (grid.n_columns != matrix.n_pixels_in_row()) {
        std::cerr << "mismatch in number of pixels with matrix\n";
        exit(-1);
      }
      if (matrix.triangular()) {
        std::cerr << "matrix is not full\n";
        exit(-1);
      }

      for (auto& element : matrix) {

        auto lor = element.lor;
        F weight = element.hits / n_emissions;
        lor_info.push_back_pixel(lor, element.pixel, weight);
      }
      lor_info.sort();
    }

    PET2D::Barrel::LMReconstruction<F, S> reconstruction(
        lor_info, cl.get<double>("sigma") / 2);
    if (verbose)
      std::cout << "created reconstruction\n";
    if (cl.exist("system"))
      reconstruction.use_system_matrix();
    if (!cl.exist("system")) {
      reconstruction.calculate_weight();
    }

    reconstruction.calculate_sensitivity();

    {
      std::ofstream out_sensitivity(output.wo_ext() + "_sensitivity" +
                                    output.ext());
      for (auto& sens : reconstruction.sensitivity())
        out_sensitivity << sens << "\n";
    }

    std::ifstream response_stream(cl.get<std::string>("response"));
    reconstruction.fscanf_responses(response_stream);
    if (verbose)
      std::cout << "read in  responses\n";

    if (cl.exist("graph")) {
      int event_num = cl.get<int>("event");
      auto graph_file_name = output.wo_ext() + ".m";
      std::ofstream graph_out(graph_file_name);

      Graphics<F> graph(graph_out);

      auto big_barrel = BarrelBuilder::make_big_barrel();
      graph.add(big_barrel);

      auto event = reconstruction.event(event_num);
      auto lor = event.lor;
      graph.add(big_barrel, lor);
      graph.add(event.p);
      for (auto it = event.first_pixel; it != event.last_pixel; ++it) {
        graph.addPixel(grid, it->pixel);
      }

      return 0;
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
                reconstruction.n_pixels * sizeof(F));
    }
#if DEBUG
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
#endif
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

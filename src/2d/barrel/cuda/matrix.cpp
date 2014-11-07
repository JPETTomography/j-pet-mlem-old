#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>
#include <exception>

#include "config.h"
#include "matrix.h"
#include "geometry.h"

#include "2d/barrel/detector_ring.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/monte_carlo.h"
#include "2d/barrel/model.h"
#include "2d/geometry/point.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

// forward declare CUDA function
bool run_matrix(Pixel<> pixel,
                int n_tof_positions,
                int number_of_threads_per_block,
                int number_of_blocks,
                int n_emissions,
                float radius,
                float h_detector,
                float w_detector,
                float pixel_size,
                unsigned int* cpu_prng_seed,
                GPU::MatrixElement* cpu_matrix,
                GPU::MatrixElement* gpu_output);

/// \brief Generated PRNG seed
/// \todo TODO: Document why are we doing that
static void gen_prng_seed(unsigned int* cpu_prng_seed,
                          int number_of_blocks,
                          int number_of_threads_per_block) {
  std::default_random_engine gen;
  std::uniform_int_distribution<unsigned int> dis(1024, 1000000);
  gen.seed(345555);
  for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block; ++i) {
    // cpu_prng_seed[i] = 53445 + i; //dis(gen);
    cpu_prng_seed[i] = dis(gen);
  }
}

GPU::OutputMatrix run_matrix(cmdline::parser& cl) {

  auto pixels_in_row = cl.get<int>("n-pixels");
  auto n_detectors = cl.get<int>("n-detectors");
  auto n_emissions = cl.get<int>("n-emissions");
  auto radius = cl.get<double>("radius");
  auto s_pixel = cl.get<double>("s-pixel");
  auto w_detector = cl.get<double>("w-detector");
  auto h_detector = cl.get<double>("h-detector");

  int* emission_adr = &cl.get<int>("n-emissions");

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,

  auto number_of_blocks = cl.get<int>("cuda-blocks");
  auto number_of_threads_per_block = cl.get<int>("cuda-threads");

  int iteration_per_thread =
      floor(n_emissions / (number_of_blocks * number_of_threads_per_block));

  printf("Gpu grid config:\n");
  printf("Number of blocks:= %d\n", number_of_blocks);
  printf("Number of threads per block:= %d\n", number_of_threads_per_block);
  printf("Number of emissions per thread:= %d\n", iteration_per_thread);

  *emission_adr =
      iteration_per_thread * number_of_blocks * number_of_threads_per_block;

  // automatic pixel size
  if (!cl.exist("radius")) {
    if (!cl.exist("s-pixel")) {
      radius = M_SQRT2;  // exact result
    } else {
      radius = s_pixel * pixels_in_row / M_SQRT2;
    }
    std::cerr << "--radius=" << radius << std::endl;
  }

  // automatic radius
  if (!cl.exist("s-pixel")) {
    if (!cl.exist("radius")) {
      s_pixel = 2. / pixels_in_row;  // exact result
    } else {
      s_pixel = M_SQRT2 * radius / pixels_in_row;
    }
    std::cerr << "--s-pixel=" << s_pixel << std::endl;
  }

  // automatic detector size
  if (!cl.exist("w-detector")) {
    w_detector = 2 * M_PI * .9 * radius / n_detectors;
    std::cerr << "--w-detector=" << w_detector << std::endl;
  }
  if (!cl.exist("h-detector")) {
    h_detector = w_detector;
    std::cerr << "--h-detector=" << h_detector << std::endl;
  }

  std::cerr << radius << " " << s_pixel << " " << w_detector << " "
            << h_detector << std::endl;

#if NO_TOF > 0
  int n_tof_positions = 1;
#else
  int n_tof_positions =
      ((int)(ceil(2 * (2 * (radius + h_detector)) / 0.01f)) + 1) / 2 * 2;
#endif

  int triangle_pixel_size = (pixels_in_row / 2 * (pixels_in_row / 2 + 1) / 2);

  std::vector<Pixel<>> lookup_table_pixel;
  lookup_table_pixel.resize(triangle_pixel_size);

  GPU::MatrixElement matrix_element;

  for (int lor_i = 0; lor_i < LORS; ++lor_i) {
    matrix_element.hit[lor_i] = 0;
  }

#if NO_TOF > 0
  std::vector<GPU::MatrixElement> gpu_vector_output;
  gpu_vector_output.resize(1);

  std::vector<GPU::MatrixElement> cpu_matrix;
  cpu_matrix.resize(number_of_blocks);
#else
  std::vector<GPU::MatrixElement> gpu_vector_output;
  gpu_vector_output.resize(n_tof_positions);

  std::vector<GPU::MatrixElement> cpu_matrix;
  cpu_matrix.resize(n_tof_positions * number_of_blocks);
#endif

  unsigned int* cpu_prng_seed =
      new unsigned int[number_of_blocks * number_of_threads_per_block * 4];
  gen_prng_seed(cpu_prng_seed, number_of_blocks, number_of_threads_per_block);

  GPU::OutputMatrix output_matrix(
      pixels_in_row, n_detectors, *emission_adr, n_tof_positions);

  double fulltime = double();

  for (Pixel<> pixel(0, 0);
       pixel < Pixel<>::end_for_n_pixels_in_row(pixels_in_row);
       ++pixel) {

    try {
      std::fill(gpu_vector_output.begin(),  //
                gpu_vector_output.end(),
                matrix_element);
      std::fill(cpu_matrix.begin(),  //
                cpu_matrix.end(),
                matrix_element);
    } catch (std::exception e) {
      std::cout << "Fill error:" << e.what() << std::endl;
    }

    std::default_random_engine gen;
    std::uniform_int_distribution<unsigned int> dis(1024, 1000000);
    gen.seed(12344 + 10 * pixel.index());

    for (int i = 0; i < 4 * number_of_blocks * number_of_threads_per_block;
         ++i) {
      cpu_prng_seed[i] = dis(gen);
    }

    clock_t start = clock();

    run_matrix(pixel,
               n_tof_positions,
               number_of_threads_per_block,
               number_of_blocks,
               iteration_per_thread,
               radius,
               h_detector,
               w_detector,
               s_pixel,
               cpu_prng_seed,
               cpu_matrix.data(),
               gpu_vector_output.data());

    clock_t stop = clock();

    fulltime += static_cast<double>(stop - start) / CLOCKS_PER_SEC;

    GPU::OutputMatrix::LOR lut_lors[LORS];
    for (GPU::OutputMatrix::LOR lor(0, 0); lor.index() < LORS; ++lor) {
      lut_lors[lor.index()] = lor;
    }

#if NO_TOF > 0
    for (int lor_i = 0; lor_i < LORS; ++lor_i) {
      if (gpu_vector_output[0].hit[lor_i] > 0) {
        output_matrix.emplace_back(
            lut_lors[lor_i], 0, pixel, gpu_vector_output[tof_i].hit[lor_i]);
      }
    }
#else
    for (int tof_i = 0; tof_i < n_tof_positions; ++tof_i) {
      for (int lor_i = 0; lor_i < LORS; ++lor_i) {
        if (gpu_vector_output[tof_i].hit[lor_i] > 0) {
          output_matrix.emplace_back(lut_lors[lor_i],
                                     tof_i,
                                     pixel,
                                     gpu_vector_output[tof_i].hit[lor_i]);
        }
      }
    }
#endif
  }

  std::cout << fulltime << " " << fulltime / triangle_pixel_size << " "
            << fulltime / triangle_pixel_size /
                   (iteration_per_thread * number_of_blocks *
                    number_of_threads_per_block) << std::endl;

  delete[] cpu_prng_seed;

  return output_matrix;
}

}  // GPU
}  // Barrel
}  // PET2D

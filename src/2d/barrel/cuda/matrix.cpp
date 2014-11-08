#include <iostream>
#include <utility>
#include <iomanip>
#include <ctime>
#include <random>
#include <vector>
#include <exception>

#include "matrix.h"

#include "2d/barrel/detector_ring.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/monte_carlo.h"
#include "2d/barrel/model.h"
#include "2d/geometry/point.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

// forward declare CUDA entry-point
bool run_matrix(Pixel pixel,
                const DetectorRing& detector_ring,
                int n_emissions,
                int n_threads_per_block,
                int n_blocks,
                float s_pixel,
                int n_positions,
                float tof_step,
                unsigned int* prng_seed,
                int* hits);

OutputMatrix run_matrix(cmdline::parser& cl) {

  DetectorRing detector_ring(cl.get<int>("n-detectors"),
                             cl.get<double>("radius"),
                             cl.get<double>("w-detector"),
                             cl.get<double>("h-detector"),
                             cl.get<double>("d-detector"));

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,

  auto n_blocks = cl.get<int>("cuda-blocks");
  auto n_threads_per_block = cl.get<int>("cuda-threads");

  auto& n_emissions = cl.get<int>("n-emissions");
  int iteration_per_thread =
      floor(n_emissions / (n_blocks * n_threads_per_block));

  printf("Gpu grid config:\n");
  printf("Number of blocks:= %d\n", n_blocks);
  printf("Number of threads per block:= %d\n", n_threads_per_block);
  printf("Number of emissions per thread:= %d\n", iteration_per_thread);

  // Number of emissions will be rounded to block size
  n_emissions = iteration_per_thread * n_blocks * n_threads_per_block;

  auto prng_seed = new unsigned int[n_blocks * n_threads_per_block * 4];
  std::default_random_engine gen;
  std::uniform_int_distribution<unsigned int> dis(1024, 1000000);
  gen.seed(345555);
  for (int i = 0; i < 4 * n_blocks * n_threads_per_block; ++i) {
    // prng_seed[i] = 53445 + i; //dis(gen);
    prng_seed[i] = dis(gen);
  }

  auto tof_step = cl.get<double>("tof-step");
  int n_tof_positions = 1;
  double max_bias = 0;
  if (cl.exist("tof-step") && tof_step > 0) {
    max_bias = Model::max_bias();
    n_tof_positions = detector_ring.n_positions(tof_step, max_bias);
  }

  auto n_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<int>("s-pixel");
  OutputMatrix output_matrix(
      n_pixels, cl.get<int>("n-detectors"), n_emissions, n_tof_positions);

  double fulltime = double();

  for (Pixel pixel(0, 0);  // start from central pixel
       pixel < Pixel::end_for_n_pixels_in_row(n_pixels);
       ++pixel) {

    std::vector<int> pixel_hits(detector_ring.n_lors * n_tof_positions, 0);

    clock_t start = clock();

    run_matrix(pixel,
               detector_ring,
               n_emissions,
               n_blocks,
               n_threads_per_block,
               s_pixel,
               n_tof_positions,
               tof_step,
               prng_seed,
               pixel_hits.data());

    clock_t stop = clock();

    fulltime += static_cast<double>(stop - start) / CLOCKS_PER_SEC;

    for (LOR lor(0, 0); lor.index() < detector_ring.n_lors; ++lor) {
      for (int position = 0; position < n_tof_positions; ++position) {
        auto hits = pixel_hits[n_tof_positions * lor.index() + position];
        if (hits > 0) {
          output_matrix.emplace_back(lor, position, pixel, hits);
        }
      }
    }
  }

  auto total_pixels = Pixel::end_for_n_pixels_in_row(n_pixels).index();
  std::cout << fulltime << " " << fulltime / total_pixels << " "
            << fulltime / total_pixels / (iteration_per_thread * n_blocks *
                                          n_threads_per_block) << std::endl;

  delete[] prng_seed;

  return output_matrix;
}

}  // GPU
}  // Barrel
}  // PET2D

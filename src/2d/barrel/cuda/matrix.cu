#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/random.h"

#include "matrix.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

__global__ static void kernel(const Pixel pixel,
                              const Scanner* scanner_ptr,
                              int n_emissions,
                              float pixel_size,
                              int n_positions,
                              float tof_step,
                              float length_scale,
                              unsigned int* gpu_prng_seed,
                              int* pixel_hits) {
  bool tof = tof_step > 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe gen(&gpu_prng_seed[4 * tid]);
  util::random::uniform_real_distribution<float> one_dis(0, 1);
  util::random::uniform_real_distribution<float> pi_dis(0, (float)M_PI);

  __shared__ util::cuda::copy<Scanner> scanner_shared_storage;
  scanner_shared_storage = scanner_ptr;
  Scanner& scanner = *scanner_shared_storage;

  Model model(length_scale);
  auto fov_radius2 = scanner.fov_radius() * scanner.fov_radius();

  for (int i = 0; i < n_emissions; ++i) {
    auto rx = (pixel.x + one_dis(gen)) * pixel_size;
    auto ry = (pixel.y + one_dis(gen)) * pixel_size;
    auto angle = pi_dis(gen);

    // ensure we are within a triangle
    if (rx > ry)
      continue;

    // ensure we are within FOV
    if (rx * rx + ry * ry > fov_radius2)
      continue;

    Event event(rx, ry, angle);
    Scanner::Response response;
    auto hits = scanner.detect(gen, model, event, response);

    int quantized_position = 0;
    if (tof)
      quantized_position =
          Scanner::quantize_tof_position(response.dl, tof_step, n_positions);

    // do we have hit on both sides?
    if (hits >= 2) {
      auto pixel_index =
          response.lor.index() * n_positions + quantized_position;
      atomicAdd(&pixel_hits[pixel_index], 1);
    }
  }

  gen.save(&gpu_prng_seed[4 * tid]);
}

Matrix<Scanner>::Matrix(const Scanner& scanner,
                        int n_threads_per_block,
                        int n_blocks,
                        float pixel_size,
                        int n_positions,
                        float tof_step,
                        float length_scale,
                        unsigned int* prng_seed)
    : n_threads_per_block(n_threads_per_block),
      n_blocks(n_blocks),
      pixel_size(pixel_size),
      n_positions(n_positions),
      tof_step(tof_step),
      length_scale(length_scale),
      pixel_hits_count(LOR::end_for_detectors(scanner.size()).index() *
                       n_positions),
      pixel_hits_size(pixel_hits_count * sizeof(int)) {

  cudaMalloc((void**)&gpu_scanner, sizeof(Scanner));
  cudaMemcpy(gpu_scanner, &scanner, sizeof(Scanner), cudaMemcpyHostToDevice);

  cudaMalloc((void**)&gpu_pixel_hits, pixel_hits_size);

  int prng_seed_size = n_blocks * n_threads_per_block * 4 * sizeof(*prng_seed);
  cudaMalloc((void**)&gpu_prng_seed, prng_seed_size);
  cudaMemcpy(gpu_prng_seed, prng_seed, prng_seed_size, cudaMemcpyHostToDevice);
}

Matrix<Scanner>::~Matrix() {
  cudaFree(gpu_prng_seed);
  cudaFree(gpu_pixel_hits);
  cudaFree(gpu_scanner);
}

void Matrix<Scanner>::operator()(Pixel pixel,
                                 int n_emissions,
                                 int* pixel_hits) {

  cudaMemset(gpu_pixel_hits, 0, pixel_hits_size);

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define kernel kernel<<<blocks, threads>>>
#endif
  kernel(pixel,
         gpu_scanner,
         n_emissions,
         pixel_size,
         n_positions,
         tof_step,
         length_scale,
         gpu_prng_seed,
         gpu_pixel_hits);

  cudaThreadSynchronize();
  cudaMemcpy(
      pixel_hits, gpu_pixel_hits, pixel_hits_size, cudaMemcpyDeviceToHost);
}

void Matrix<Scanner>::run(Scanner& scanner,
                          int n_blocks,
                          int n_threads_per_block,
                          int n_emissions,
                          double tof_step,
                          int n_pixels,
                          double s_pixel,
                          double length_scale,
                          util::delegate<void(int, bool)> progress,
                          util::delegate<void(LOR, S, Pixel, Hit)> entry) {

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,
  auto n_threads = n_blocks * n_threads_per_block;
  auto n_thread_emissions = (n_emissions + n_threads - 1) / n_threads;
  // Number of emissions will be rounded to block size
  n_emissions = n_thread_emissions * n_blocks * n_threads_per_block;

  int n_tof_positions = 1;
  double max_bias = 0;
  if (tof_step > 0) {
    max_bias = 0;
    n_tof_positions = scanner.n_tof_positions(tof_step, max_bias);
  }

  unsigned int prng_seed[n_blocks * n_threads_per_block * 4];
  util::random::tausworthe gen;
  gen.seed(345555);
  for (int i = 0; i < 4 * n_blocks * n_threads_per_block; ++i) {
    prng_seed[i] = gen();  // 53445 + i
  }

  GPU::Matrix<Scanner> gpu_matrix(scanner,
                                  n_threads_per_block,
                                  n_blocks,
                                  s_pixel,
                                  n_tof_positions,
                                  tof_step,
                                  length_scale,
                                  prng_seed);

  const auto n_lors = LOR::end_for_detectors(scanner.size()).index();
  int pixel_hits[n_lors * n_tof_positions];

  LOR lor_map[n_lors];
  for (LOR lor(0, 0); lor < LOR::end_for_detectors(scanner.size()); ++lor) {
    lor_map[lor.index()] = lor;
  }

  auto end_pixel = Pixel::end_for_n_pixels_in_row(n_pixels / 2);

  for (Pixel pixel(0, 0); pixel < end_pixel; ++pixel) {
    progress(pixel.index(), false);

    gpu_matrix(pixel, n_thread_emissions, pixel_hits);

    for (size_t lor_index = 0; lor_index < n_lors; ++lor_index) {
      auto lor = lor_map[lor_index];
      for (int position = 0; position < n_tof_positions; ++position) {
        auto hits = pixel_hits[n_tof_positions * lor_index + position];
        if (hits > 0) {
          entry(lor, position, pixel, hits);
        }
      }
    }

    progress(pixel.index(), true);
  }
}

}  // GPU
}  // Barrel
}  // PET2D

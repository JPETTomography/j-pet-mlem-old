#include <cuda_runtime.h>
#include <stdio.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/random.h"

#include "matrix.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
namespace Matrix {

__global__ static void kernel(const Pixel pixel,
                              const Scanner* scanner_ptr,
                              int n_emissions,
                              float s_pixel,
                              int n_tof_positions,
                              float tof_step,
                              float length_scale,
                              unsigned int* gpu_rng_seed,
                              int* pixel_hits) {
  bool tof = tof_step > 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe rng(&gpu_rng_seed[4 * tid]);
  util::random::uniform_real_distribution<float> one_dis(0, 1);
  util::random::uniform_real_distribution<float> pi_dis(0, (float)M_PI);

  __shared__ util::cuda::copy<Scanner> scanner_shared_storage;
  scanner_shared_storage = scanner_ptr;
  Scanner& scanner = *scanner_shared_storage;

  Model model(length_scale);
  auto fov_radius2 = scanner.fov_radius() * scanner.fov_radius();

  for (int i = 0; i < n_emissions; ++i) {
    auto rx = (pixel.x + one_dis(rng)) * s_pixel;
    auto ry = (pixel.y + one_dis(rng)) * s_pixel;
    auto angle = pi_dis(rng);

    // ensure we are within a triangle
    if (rx > ry)
      continue;

    // ensure we are within FOV
    if (rx * rx + ry * ry > fov_radius2)
      continue;

    Event event(rx, ry, angle);
    Scanner::Response response;
    auto hits = scanner.detect(rng, model, event, response);

    int quantized_position = 0;
    if (tof)
      quantized_position = Scanner::quantize_tof_position(
          response.dl, tof_step, n_tof_positions);

    // do we have hit on both sides?
    if (hits >= 2) {
      auto pixel_index =
          response.lor.index() * n_tof_positions + quantized_position;
      atomicAdd(&pixel_hits[pixel_index], 1);
    }
  }

  rng.save(&gpu_rng_seed[4 * tid]);
}

template <>
void run<Scanner>(Scanner& scanner,
                  util::random::tausworthe& rng,
                  int n_blocks,
                  int n_threads_per_block,
                  int n_emissions,
                  double tof_step,
                  int n_tof_positions,
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
  n_emissions = n_thread_emissions * n_threads;

  // (1) setup LORs:
  const auto end_lor = LOR::end_for_detectors(scanner.size());
  const auto n_lors = end_lor.index();
  LOR lor_map[n_lors];
  for (LOR lor(0, 0); lor < end_lor; ++lor) {
    lor_map[lor.index()] = lor;
  }

  // (2) copy scanner to GPU:
  Scanner* gpu_scanner;
  cudaMalloc((void**)&gpu_scanner, sizeof(Scanner));
  cudaMemcpy(gpu_scanner, &scanner, sizeof(Scanner), cudaMemcpyHostToDevice);

  // (3) create empty pixel hits for GPU:
  const int pixel_hits_len = n_lors * n_tof_positions;
  const int pixel_hits_size = pixel_hits_len * sizeof(int);
  int* gpu_pixel_hits;
  int pixel_hits[pixel_hits_len];
  cudaMalloc((void**)&gpu_pixel_hits, pixel_hits_size);

  // (4) initalize RNG for all GPU threads:
  int rng_seed_len = n_threads * 4;
  unsigned int rng_seed[rng_seed_len];
  int rng_seed_size = rng_seed_len * sizeof(*rng_seed);
  for (int i = 0; i < rng_seed_len; ++i) {
    rng_seed[i] = rng();
  }
  unsigned int* gpu_rng_seed;
  cudaMalloc((void**)&gpu_rng_seed, rng_seed_size);
  cudaMemcpy(gpu_rng_seed, rng_seed, rng_seed_size, cudaMemcpyHostToDevice);

  auto end_pixel = Pixel::end_for_n_pixels_in_row(n_pixels / 2);
  for (Pixel pixel(0, 0); pixel < end_pixel; ++pixel) {
    progress(pixel.index(), false);

    // (5) clear pixel hits on GPU:
    cudaMemset(gpu_pixel_hits, 0, pixel_hits_size);

#if __CUDACC__
    dim3 blocks(n_blocks);
    dim3 threads(n_threads_per_block);
#define kernel kernel<<<blocks, threads>>>
#endif
    // (6) run kernel:
    kernel(pixel,
           gpu_scanner,
           n_emissions,
           s_pixel,
           n_tof_positions,
           tof_step,
           length_scale,
           gpu_rng_seed,
           gpu_pixel_hits);

    // (7) copy back the data:
    cudaThreadSynchronize();
    cudaMemcpy(
        pixel_hits, gpu_pixel_hits, pixel_hits_size, cudaMemcpyDeviceToHost);

    // (8) inform callee about data:
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

  // (9) free the GPU data:
  cudaFree(gpu_rng_seed);
  cudaFree(gpu_pixel_hits);
  cudaFree(gpu_scanner);
}

}  // Matrix
}  // GPU
}  // Barrel
}  // PET2D

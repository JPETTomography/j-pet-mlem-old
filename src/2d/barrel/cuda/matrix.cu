#include <cuda_runtime.h>

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
                              int n_thread_emissions,
                              F s_pixel,
                              int n_tof_positions,
                              F tof_step,
                              F length_scale,
                              util::random::tausworthe::state_type* rng_state,
                              int* pixel_hits) {
  bool tof = tof_step > 0;
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  util::random::tausworthe rng(rng_state[tid]);
  util::random::uniform_real_distribution<F> one_dis(0, 1);
  util::random::uniform_real_distribution<F> pi_dis(0, (F)M_PI);

  __shared__ util::cuda::copy<Scanner> scanner_shared_storage;
  scanner_shared_storage = scanner_ptr;
  Scanner& scanner = *scanner_shared_storage;

  Model model(length_scale);

  for (int i = 0; i < n_thread_emissions; ++i) {
    auto rx = (pixel.x + one_dis(rng)) * s_pixel;
    auto ry = (pixel.y + one_dis(rng)) * s_pixel;

    // ensure we are within a triangle, so we got only half hits on diagonal
    if (rx > ry)
      continue;

    auto angle = pi_dis(rng);

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

  rng.save(rng_state[tid]);
}

template <>
void run<Scanner>(
    Scanner& scanner,
    util::random::tausworthe& rng,
    int n_emissions,
    double tof_step,
    int n_tof_positions,
    int n_pixels,
    double s_pixel,
    double length_scale,
    util::delegate<void(int completed, bool finished)> progress,
    util::delegate<void(LOR lor, S position, Pixel pixel, Hit hits)> entry,
    int device,
    int n_blocks,
    int n_threads_per_block,
    util::delegate<void(const char* device_name, int n_emissions)> info) {

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define kernel kernel<<<blocks, threads>>>
#endif

  // GTX 770 - 8 SMX * 192 cores = 1536 cores -
  // each SMX can use 8 active blocks,
  auto n_threads = n_blocks * n_threads_per_block;
  auto n_thread_emissions = (n_emissions + n_threads - 1) / n_threads;
  // Number of emissions will be rounded to block size
  n_emissions = n_thread_emissions * n_threads;

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  info(prop.name, n_emissions);

  const auto end_lor = LOR::end_for_detectors(scanner.size());
  const auto n_lors = end_lor.index();
  LOR lor_map[n_lors];
  for (LOR lor(0, 0); lor < end_lor; ++lor) {
    lor_map[lor.index()] = lor;
  }

  util::cuda::on_device<Scanner> scanner_on_device(scanner);
  util::cuda::memory<int> pixel_hits(n_lors * n_tof_positions);
  util::cuda::memory<util::random::tausworthe::state_type> rng_state(n_threads);

  for (size_t i = 0; i < rng_state.size; ++i) {
    util::random::tausworthe thread_rng(rng);
    thread_rng.save(rng_state[i]);
  }
  rng_state.copy_to_device();

  const auto pixel_fov_radius = scanner.fov_radius() / s_pixel;
  const int pixel_fov_radius2 = pixel_fov_radius * pixel_fov_radius;

  auto end_pixel = Pixel::end_for_n_pixels_in_row(n_pixels / 2);
  for (Pixel pixel(0, 0); pixel < end_pixel; ++pixel) {
    // ensure we are within FOV
    if (pixel.distance_from_origin2() > pixel_fov_radius2)
      continue;

    progress(pixel.index(), false);

    pixel_hits.zero_on_device();

    kernel(pixel,
           scanner_on_device,
           n_thread_emissions,
           s_pixel,
           n_tof_positions,
           tof_step,
           length_scale,
           rng_state,
           pixel_hits);
    cudaThreadSynchronize();

    pixel_hits.copy_from_device();

    for (int lor_index = 0; lor_index < n_lors; ++lor_index) {
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
  progress(end_pixel.index() - 1, true);
}

}  // Matrix
}  // GPU
}  // Barrel
}  // PET2D

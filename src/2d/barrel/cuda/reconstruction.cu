#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "common/cuda/kernels.h"

#include "reconstruction.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
namespace Reconstruction {

texture<F, 2, cudaReadModeElementType> tex_rho;

// foreach p: count y[p] and store it in output_rho[p]
__global__ static void reconstruction_1(const Pixel* pixels,
                                        const F* pixel_weights,
                                        const size_t* lor_pixel_info_begin,
                                        const size_t* lor_pixel_info_end,
                                        const Mean* means,
                                        const int n_means,
                                        F* output_rho,
                                        const int width) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_means + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int mean_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (mean_index >= n_means) {
      break;
    }

    auto mean = means[mean_index];
    auto lor_index = mean.lor.index();
    auto pixel_info_begin = lor_pixel_info_begin[lor_index];
    auto pixel_info_end = lor_pixel_info_end[lor_index];

    // count u for current lor
    F u = 0;
    for (auto info_index = pixel_info_begin; info_index < pixel_info_end;
         ++info_index) {
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      u += tex2D(tex_rho, pixel.x, pixel.y) * pixel_weight;
    }
    F phi = mean.mean / u;
    for (auto info_index = pixel_info_begin; info_index < pixel_info_end;
         ++info_index) {
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      atomicAdd(&output_rho[pixel.index(width)], phi * pixel_weight);
    }
  }
}

// foreach p: count output_rho[p] *= rho[p]
__global__ static void reconstruction_2(F* output_rho,
                                        const F* scale,
                                        const int width,
                                        const int height) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_pixels = width * height;
  const auto n_pixel_chunks = (n_pixels + n_threads - 1) / n_threads;

  for (int chunk = 0; chunk < n_pixel_chunks; ++chunk) {
    int pixel_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (pixel_index >= n_pixels) {
      break;
    }
    Pixel pixel(pixel_index % width, pixel_index / width);
    // there is no collision there, so we don't need atomics
    output_rho[pixel_index] *=
        tex2D(tex_rho, pixel.x, pixel.y) * scale[pixel_index];
  }
}

void run(const Geometry& geometry,
         const Mean* means,
         int n_means,
         int width,
         int height,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, const Output& output)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> info) {

#if __CUDACC__
  dim3 blocks(n_blocks);
  dim3 threads(n_threads_per_block);
#define reduce_to_sensitivity reduce_to_sensitivity<<<blocks, threads>>>
#define invert invert<<<blocks, threads>>>
#define reconstruction_1 reconstruction_1<<<blocks, threads>>>
#define reconstruction_2 reconstruction_2<<<blocks, threads>>>
#else
  (void)n_blocks, n_threads_per_block;  // mark used
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  info(prop.name);

  util::cuda::on_device<Pixel> device_pixels(geometry.pixels,
                                             geometry.n_pixel_infos);
  util::cuda::on_device<F> device_pixel_weights(geometry.pixel_weights,
                                                geometry.n_pixel_infos);
  util::cuda::on_device<size_t> device_lor_pixel_info_begin(
      geometry.lor_pixel_info_begin, geometry.n_lors);
  util::cuda::on_device<size_t> device_lor_pixel_info_end(
      geometry.lor_pixel_info_end, geometry.n_lors);
  util::cuda::on_device<Mean> device_means(means, n_means);

  util::cuda::texture2D<F> rho(tex_rho, width, height);
  util::cuda::memory<F> output_rho((size_t)width * height);
  Output rho_output(width, height, output_rho.host_ptr);

  for (auto& v : output_rho) {
    v = 1;
  }
  output_rho.copy_to_device();

  util::cuda::on_device<F> scale((size_t)width * height);
  scale.zero_on_device();

  Common::GPU::reduce_to_sensitivity(device_pixels,
                                     device_pixel_weights,
                                     geometry.n_pixel_infos,
                                     scale,
                                     width);
  cudaThreadSynchronize();

  Common::GPU::invert(scale, width * height);
  cudaThreadSynchronize();

  for (int block = 0; block < n_iteration_blocks; ++block) {
    for (int i = 0; i < n_iterations_in_block; ++i) {
      int iteration = block * n_iterations_in_block + i;
      progress(iteration, false);

      rho = output_rho;
      output_rho.zero_on_device();

      reconstruction_1(device_pixels,
                       device_pixel_weights,
                       device_lor_pixel_info_begin,
                       device_lor_pixel_info_end,
                       device_means,
                       n_means,
                       output_rho,
                       width);
      cudaThreadSynchronize();

      reconstruction_2(output_rho, scale, width, height);
      cudaThreadSynchronize();

      progress(iteration, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!block && i < n_iterations_in_block - 1 &&
          (i < 5 || i == 9 || i == 14 || i == 19 || i == 29 || i == 49 ||
           i == 99)) {
        output_rho.copy_from_device();
        output(i + 1, rho_output);
      }
    }

    output_rho.copy_from_device();
    int iteration = (block + 1) * n_iterations_in_block;
    output(iteration, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D

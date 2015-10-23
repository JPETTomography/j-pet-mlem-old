#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "common/cuda/kernels.h"

#include "lm_reconstruction.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
namespace LMReconstruction {

texture<float, 2, cudaReadModeElementType> tex_rho;

template <typename FType, typename SType> struct Kernel {
  using F = FType;
  using S = SType;

  __device__ Kernel(const F* scale, F sigma, S width)
      : scale(scale),
        sigma(sigma),
        gauss_norm_dl(1 / (sigma * compat::sqrt(2 * M_PI))),
        inv_sigma2_dl(1 / (2 * sigma * sigma)),
        width(width) {}

  __device__ F l(const Event& event, const PixelInfo& pixel_info) {
    auto diff_t = pixel_info.t - event.t;
    return gauss_norm_dl * compat::exp(-diff_t * diff_t * inv_sigma2_dl);
  }

  __device__ F operator()(const Event& event, const PixelInfo& pixel_info) {
    const auto pixel = pixel_info.pixel;
    const auto pixel_index = pixel.y * width + pixel.x;
    const auto kernel_z = pixel_info.weight * scale[pixel_index];
    return l(event, pixel_info) * kernel_z * tex2D(tex_rho, pixel.x, pixel.y);
  }

  const F* scale;
  const F sigma;
  const F gauss_norm_dl;
  const F inv_sigma2_dl;
  const S width;
};

__global__ static void reconstruction(const PixelInfo* pixel_infos,
                                      const Event* events,
                                      const int n_events,
                                      float* output_rho,
                                      const float* scale,
                                      const float sigma,
                                      const int width) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_events + n_threads - 1) / n_threads;

  Kernel<float, int> kernel(scale, sigma, width);

  // --- event loop ----------------------------------------------------------
  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int event_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (event_index >= n_events) {
      break;
    }

    const auto event = events[event_index];
    F denominator = 0;

    // -- voxel loop - denominator -------------------------------------------
    for (auto info_index = event.first_pixel_info_index;
         info_index < event.last_pixel_info_index;
         ++info_index) {
      const auto& pixel_info = pixel_infos[info_index];
      denominator += kernel(event, pixel_info);
    }  // voxel loop - denominator

    if (denominator == 0)
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto info_index = event.first_pixel_info_index;
         info_index < event.last_pixel_info_index;
         ++info_index) {
      const auto& pixel_info = pixel_infos[info_index];
      const auto pixel_index = pixel_info.pixel.y * width + pixel_info.pixel.x;
      atomicAdd(&output_rho[pixel_index],
                kernel(event, pixel_info) * inv_denominator);
    }  // voxel loop
  }    // event loop
}

__global__ static void add_offsets(Event* events,
                                   const int n_events,
                                   const size_t* lor_pixel_info_start) {
  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_events + n_threads - 1) / n_threads;
  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int event_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (event_index >= n_events) {
      break;
    }
    auto& event = events[event_index];
    const auto pixel_info_start = lor_pixel_info_start[event.lor.index()];
    event.first_pixel_info_index += pixel_info_start;
    event.last_pixel_info_index += pixel_info_start;
  }
}

void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         float sigma,
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
#define sensitivity sensitivity<<<blocks, threads>>>
#define invert invert<<<blocks, threads>>>
#define reconstruction reconstruction<<<blocks, threads>>>
#define add_offsets add_offsets<<<blocks, threads>>>
#else
  (void)n_blocks, n_threads_per_block;  // mark used
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  info(prop.name);

  util::cuda::on_device<PixelInfo> device_pixel_infos(geometry.pixel_infos,
                                                      geometry.n_pixel_infos);
  util::cuda::on_device<size_t> device_lor_pixel_info_start(
      geometry.lor_pixel_info_start, geometry.n_lors);
  util::cuda::on_device<Event> device_events(events, n_events);

  add_offsets(device_events, n_events, device_lor_pixel_info_start);

  util::cuda::memory2D<F> rho(tex_rho, width, height);
  util::cuda::on_device<F> output_rho((size_t)width * height);

  for (auto& v : rho) {
    v = 1;
  }
  rho.copy_to_device();

  util::cuda::on_device<F> scale((size_t)width * height);
  scale.zero_on_device();

  Common::GPU::sensitivity(
      device_pixel_infos, geometry.n_pixel_infos, scale, width);
  cudaThreadSynchronize();

  Common::GPU::invert(scale, width * height);
  cudaThreadSynchronize();

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      output_rho.zero_on_device();

      reconstruction(device_pixel_infos,
                     device_events,
                     n_events,
                     output_rho,
                     scale,
                     sigma,
                     width);
      cudaThreadSynchronize();

      rho = output_rho;
      progress(ib * n_iterations_in_block + it, true);
    }

    rho.copy_from_device();
    Output rho_output(width, height, rho.host_ptr);
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D

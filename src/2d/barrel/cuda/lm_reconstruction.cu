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

texture<F, 2, cudaReadModeElementType> tex_rho;

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
                                      F* output_rho,
                                      const F* scale,
                                      const F sigma,
                                      const int width) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_events + n_threads - 1) / n_threads;

  Kernel<F, int> kernel(scale, sigma, width);

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
    for (auto info_index = event.pixel_info_begin;
         info_index < event.pixel_info_end;
         ++info_index) {
      const auto& pixel_info = pixel_infos[info_index];
      denominator += kernel(event, pixel_info);
    }  // voxel loop - denominator

    if (denominator == 0)
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto info_index = event.pixel_info_begin;
         info_index < event.pixel_info_end;
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
                                   const size_t* lor_pixel_info_begin) {
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
    const auto pixel_info_begin = lor_pixel_info_begin[event.lor.index()];
    event.pixel_info_begin += pixel_info_begin;
    event.pixel_info_end += pixel_info_begin;
  }
}

void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         F sigma,
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
  util::cuda::on_device<size_t> device_lor_pixel_info_begin(
      geometry.lor_pixel_info_begin, geometry.n_lors);
  util::cuda::on_device<Event> device_events(events, n_events);

  add_offsets(device_events, n_events, device_lor_pixel_info_begin);

  util::cuda::texture2D<F> rho(tex_rho, width, height);
  util::cuda::memory<F> output_rho((size_t)width * height);

  for (auto& v : output_rho) {
    v = 1;
  }
  output_rho.copy_to_device();

  util::cuda::on_device<F> scale((size_t)width * height);
  scale.zero_on_device();

  Common::GPU::reduce_to_sensitivity(
      device_pixel_infos, geometry.n_pixel_infos, scale, width);
  cudaThreadSynchronize();

  Common::GPU::invert(scale, width * height);
  cudaThreadSynchronize();

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      rho = output_rho;
      output_rho.zero_on_device();

      reconstruction(device_pixel_infos,
                     device_events,
                     n_events,
                     output_rho,
                     scale,
                     sigma,
                     width);
      cudaThreadSynchronize();

      progress(ib * n_iterations_in_block + it, true);
    }

    output_rho.copy_from_device();
    Output rho_output(width, height, output_rho.host_ptr);
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D

#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "reconstruction.h"

#define USE_TEXTURE 1  // using textures is faster when using bigger rhos

#if USE_WARP_GRANULARITY
#if USE_VOXEL_GRANULARITY  // voxel (densier) granularity
#include "reconstruction/warp_voxel_granularity.cuh"
#else            // pixel granularity (faster)
#if USE_TEXTURE  // use texture and 3D arrays for rho and 2D sensitivity lookup
#include "reconstruction/warp_granularity.cuh"
#else  // use linear memory
#include "reconstruction/warp_granularity_no_tex.cuh"
#endif
#endif
#elif USE_THREAD_GRANULARITY
#include "reconstruction/thread_granularity.cuh"
#endif

#include "common/cuda/kernels.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {
namespace Reconstruction {

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
         const Sensitivity& sensitivity,
         const Event* events,
         int n_events,
         F sigma_z,
         F sigma_dl,
         const Grid& grid,
         const F barrel_length,
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

  util::cuda::on_device<LineSegment> device_lor_line_segments(
      geometry.lor_line_segments, geometry.n_lors);
  util::cuda::on_device<PixelInfo> device_pixel_infos(geometry.pixel_infos,
                                                      geometry.n_pixel_infos);
  util::cuda::on_device<size_t> device_lor_pixel_info_begin(
      geometry.lor_pixel_info_begin, geometry.n_lors);
  util::cuda::on_device<Event> device_events(events, n_events);

  add_offsets(device_events, n_events, device_lor_pixel_info_begin);

#if USE_TEXTURE
  util::cuda::texture3D<F> rho(tex_rho,
                               grid.pixel_grid.n_columns,
                               grid.pixel_grid.n_rows,
                               grid.n_planes);
#else
  util::cuda::on_device<F> rho((size_t)grid.n_voxels);
#endif
  util::cuda::memory<F> output_rho((size_t)grid.n_voxels);
  Output rho_output(grid.pixel_grid.n_columns,
                    grid.pixel_grid.n_rows,
                    grid.n_planes,
                    output_rho.host_ptr);

  for (auto& v : output_rho) {
    v = 1;
  }
  output_rho.copy_to_device();

#if USE_TEXTURE
  util::cuda::texture2D<F> device_sensitivity(tex_sensitivity,
                                              (size_t)sensitivity.width,
                                              (size_t)sensitivity.height,
                                              sensitivity.data);
#else
  util::cuda::on_device<F> device_sensitivity((size_t)grid.pixel_grid.n_pixels);
#endif
  (void)device_sensitivity;  // device sensitivity is used via tex_sensitivity

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      rho = output_rho;
      output_rho.zero_on_device();

      reconstruction(device_lor_line_segments,
                     device_pixel_infos,
                     device_events,
                     n_events,
                     output_rho,
#if !USE_TEXTURE
                     rho,
                     device_sensitivity,
#endif
                     sigma_z,
                     sigma_dl,
                     grid,
                     barrel_length);
      cudaThreadSynchronize();

      progress(ib * n_iterations_in_block + it, true);

      // always output first 5 iterations, and at 10, 15, 20, 30, 50, 100
      if (!ib && it < n_iterations_in_block - 1 &&
          (it < 5 || it == 9 || it == 14 || it == 19 || it == 29 || it == 49 ||
           it == 99)) {
        output_rho.copy_from_device();
        output(it + 1, rho_output);
      }
    }

    output_rho.copy_from_device();
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

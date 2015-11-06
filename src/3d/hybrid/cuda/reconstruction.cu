#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "common/cuda/kernels.h"

#include "reconstruction.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {
namespace Reconstruction {

texture<float, 3, cudaReadModeElementType> tex_rho;

__global__ static void reconstruction(const PixelInfo* pixel_infos,
                                      const Event* events,
                                      const int n_events,
                                      float* output_rho,
                                      const float* scale,
                                      const float sigma_z,
                                      const float sigma_dl,
                                      const int width,
                                      const int height) {

  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_threads = gridDim.x * blockDim.x;
  const auto n_chunks = (n_events + n_threads - 1) / n_threads;

  // --- event loop ----------------------------------------------------------
  for (int chunk = 0; chunk < n_chunks; ++chunk) {
    int event_index = chunk * n_threads + tid;
    // check if we are still on the list
    if (event_index >= n_events) {
      break;
    }

    const auto event = events[event_index];
    (void)event;  // FIXME: implement me!

  }  // event loop
}

void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         float sigma_z,
         float sigma_dl,
         int width,
         int height,
         int depth,
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

  util::cuda::memory3D<F> rho(tex_rho, width, height, depth);
  util::cuda::on_device<F> output_rho((size_t)width * height * depth);

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
                     sigma_z,
                     sigma_dl,
                     width,
                     height);
      cudaThreadSynchronize();

      rho = output_rho;
      progress(ib * n_iterations_in_block + it, true);
    }

    rho.copy_from_device();
    Output rho_output(width, height, depth, rho.host_ptr);
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

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
#define USE_3D_SENSITIVITY 1
#include "reconstruction/warp_granularity.cuh"
#define USE_MULTI_PLANE_WEIGHTS 1
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

void run(const Geometry& geometry,
         const Output& sensitivity,
         const Event* events,
         int n_events,
         F sigma_z,
         F sigma_dl,
         const Grid& grid,
         const F barrel_length,
         int n_iteration_blocks,
         int n_iterations_in_block,
         const Output& rho_output,
         const int start_iteration,
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
#else
  (void)n_blocks, n_threads_per_block;  // mark used
#endif

  cudaSetDevice(device);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, device);
  info(prop.name);

  util::cuda::on_device<LineSegment> device_lor_line_segments(
      geometry.lor_line_segments, geometry.n_lors);
  util::cuda::on_device<Pixel> device_pixels(geometry.pixels,
                                             geometry.n_pixel_infos);
  util::cuda::on_device<F> device_pixel_weights(
      geometry.pixel_weights, geometry.n_pixel_infos * geometry.n_planes_half);
  util::cuda::on_device<size_t> device_lor_pixel_info_begin(
      geometry.lor_pixel_info_begin, geometry.n_lors);
  util::cuda::on_device<Event> device_events(events, n_events);

#if USE_TEXTURE
  util::cuda::texture3D<F> rho(tex_rho,
                               grid.pixel_grid.n_columns,
                               grid.pixel_grid.n_rows,
                               grid.n_planes);
#else
  util::cuda::on_device<F> rho((size_t)grid.n_voxels);
#endif
  util::cuda::memory<F> output_rho((size_t)grid.n_voxels, rho_output.data);
  output_rho.copy_to_device();

#if USE_TEXTURE
  util::cuda::texture2D<F> device_sensitivity(tex_sensitivity,
                                              (size_t)sensitivity.width,
                                              (size_t)sensitivity.height,
                                              sensitivity.data);
  (void)device_sensitivity;  // device sensitivity is used via tex_sensitivity
#if USE_3D_SENSITIVITY || USE_MULTI_PLANE_WEIGHTS
  util::cuda::texture3D<F> device_3d_sensitivity(tex_3d_sensitivity,
                                                 (size_t)sensitivity.width,
                                                 (size_t)sensitivity.height,
                                                 (size_t)sensitivity.depth,
                                                 sensitivity.data);
  (void)device_3d_sensitivity;  // device 3D sensitivity is
                                // used via tex_3d_sensitivity
#endif
#else
  util::cuda::on_device<F> device_sensitivity((size_t)grid.pixel_grid.n_pixels);
#endif

  for (int block = 0; block < n_iteration_blocks; ++block) {
    for (int i = 0; i < n_iterations_in_block; ++i) {
      int iteration = block * n_iterations_in_block + i;
      if (iteration < start_iteration)
        continue;
      progress(iteration, false);

      rho = output_rho;
      output_rho.zero_on_device();

#if USE_MULTI_PLANE_WEIGHTS
      if (geometry.n_planes_half > 1)
        Multiplane::reconstruction(device_lor_line_segments,
                                   device_pixels,
                                   device_pixel_weights,
                                   geometry.n_pixel_infos,
                                   geometry.n_planes_half,
                                   device_events,
                                   n_events,
                                   output_rho,
                                   sigma_z,
                                   sigma_dl,
                                   grid,
                                   barrel_length);
      else
#endif
#if USE_3D_SENSITIVITY
          if (sensitivity.depth > 1)
        Sensitivity3D::reconstruction(device_lor_line_segments,
                                      device_pixels,
                                      device_pixel_weights,
                                      sensitivity.depth,
                                      device_events,
                                      n_events,
                                      output_rho,
                                      sigma_z,
                                      sigma_dl,
                                      grid,
                                      barrel_length);
      else
#endif
        reconstruction(device_lor_line_segments,
                       device_pixels,
                       device_pixel_weights,
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
      cudaPeekAtLastError();  // ensure kernel was run successfully
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

    int iteration = (block + 1) * n_iterations_in_block;
    if (iteration <= start_iteration)
      continue;
    output_rho.copy_from_device();
    output(iteration, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

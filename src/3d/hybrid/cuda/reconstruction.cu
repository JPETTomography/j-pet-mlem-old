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

texture<F, 3, cudaReadModeElementType> tex_rho;

__global__ static void reconstruction(const LineSegment* lor_line_segments,
                                      const PixelInfo* pixel_infos,
                                      const Event* events,
                                      const int n_events,
                                      F* output_rho,
                                      const F* scale,
                                      const F sigma_z,
                                      const F sigma_dl,
                                      const Grid grid) {

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
    const auto lor = event.lor;
    const auto lor_index = lor.index();
    const auto segment = lor_line_segments[lor_index];
    const auto R = segment.length / 2;
    F denominator = 0;

    // -- voxel loop - denominator -------------------------------------------
    for (auto info_index = event.first_pixel_info_index;
         info_index < event.last_pixel_info_index;
         ++info_index) {
      const auto& pixel_info = pixel_infos[info_index];
      auto pixel = pixel_info.pixel;
      auto pixel_index = grid.pixel_grid.index(pixel);
      auto center =
          grid.pixel_grid.center_at(pixel);  // FIXME: wrong!! temporary
      auto up = segment.projection_relative_middle(center);

      for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
      }
    }  // voxel loop - denominator

    if (denominator == 0)
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto info_index = event.first_pixel_info_index;
         info_index < event.last_pixel_info_index;
         ++info_index) {
      const auto& pixel_info = pixel_infos[info_index];
      auto pixel = pixel_info.pixel;
    }  // voxel loop
  }    // event loop
}

void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         F sigma_z,
         F sigma_dl,
         const Grid& grid,
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

  util::cuda::on_device<LineSegment> device_lor_line_segments(
      geometry.lor_line_segments, geometry.n_lors);
  util::cuda::on_device<PixelInfo> device_pixel_infos(geometry.pixel_infos,
                                                      geometry.n_pixel_infos);
  util::cuda::on_device<size_t> device_lor_pixel_info_start(
      geometry.lor_pixel_info_start, geometry.n_lors);
  util::cuda::on_device<Event> device_events(events, n_events);

  util::cuda::memory3D<F> rho(tex_rho,
                              grid.pixel_grid.n_columns,
                              grid.pixel_grid.n_rows,
                              grid.n_planes);
  util::cuda::on_device<F> output_rho((size_t)grid.n_voxels);

  for (auto& v : rho) {
    v = 1;
  }
  rho.copy_to_device();

  util::cuda::on_device<F> scale((size_t)grid.pixel_grid.n_pixels);
  scale.zero_on_device();

  Common::GPU::sensitivity(device_pixel_infos,
                           geometry.n_pixel_infos,
                           scale,
                           grid.pixel_grid.n_columns);
  cudaThreadSynchronize();

  Common::GPU::invert(scale, grid.pixel_grid.n_pixels);
  cudaThreadSynchronize();

  for (int ib = 0; ib < n_iteration_blocks; ++ib) {
    for (int it = 0; it < n_iterations_in_block; ++it) {
      progress(ib * n_iterations_in_block + it, false);

      output_rho.zero_on_device();

      reconstruction(device_lor_line_segments,
                     device_pixel_infos,
                     device_events,
                     n_events,
                     output_rho,
                     scale,
                     sigma_z,
                     sigma_dl,
                     grid);
      cudaThreadSynchronize();

      rho = output_rho;
      progress(ib * n_iterations_in_block + it, true);
    }

    rho.copy_from_device();
    Output rho_output(grid.pixel_grid.n_columns,
                      grid.pixel_grid.n_rows,
                      grid.n_planes,
                      rho.host_ptr);
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

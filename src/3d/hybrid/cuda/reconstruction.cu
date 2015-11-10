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
texture<F, 2, cudaReadModeElementType> tex_sensitivity;

__global__ static void reconstruction(const LineSegment* lor_line_segments,
                                      const PixelInfo* pixel_infos,
                                      const Event* events,
                                      const int n_events,
                                      F* output_rho,
                                      const F sigma_z,
                                      const F sigma_dl,
                                      const Grid grid) {

  const Kernel kernel(sigma_z, sigma_dl);
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
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);

      for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        auto z = grid.center_z_at(voxel);
        auto diff = Point2D(up, z) - Point2D(event.up, event.right);
        auto kernel2d =
            kernel(event.up, event.tan, event.sec, R, Vector2D(diff.y, diff.x));
        auto kernel_t = pixel_info.weight;
        auto weight = kernel2d * kernel_t *  // hybrid of 2D x-y & y-z
                      tex3D(tex_rho, voxel.x, voxel.y, voxel.z);
        // end of kernel calculation
        denominator += weight * tex2D(tex_sensitivity, voxel.x, voxel.y);
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
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);

      for (int iz = event.first_plane; iz < event.last_plane; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        auto z = grid.center_z_at(voxel);
        auto diff = Point2D(up, z) - Point2D(event.up, event.right);
        auto kernel2d =
            kernel(event.up, event.tan, event.sec, R, Vector2D(diff.y, diff.x));
        auto kernel_t = pixel_info.weight;
        auto weight = kernel2d * kernel_t *  // hybrid of 2D x-y & y-z
                      tex3D(tex_rho, voxel.x, voxel.y, voxel.z);
        // end of kernel calculation

        int voxel_index = grid.index(voxel);
        atomicAdd(&output_rho[voxel_index], weight * inv_denominator);
      }
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

  util::cuda::on_device3D<F> rho(tex_rho,
                                 grid.pixel_grid.n_columns,
                                 grid.pixel_grid.n_rows,
                                 grid.n_planes);
  util::cuda::memory<F> output_rho((size_t)grid.n_voxels);

  for (auto& v : output_rho) {
    v = 1;
  }
  output_rho.copy_to_device();

  util::cuda::memory2D<F> sensitivity2d(
      tex_sensitivity, grid.pixel_grid.n_columns, grid.pixel_grid.n_rows);
  sensitivity2d.zero_on_device();

  Common::GPU::sensitivity(device_pixel_infos,
                           geometry.n_pixel_infos,
                           sensitivity2d,
                           grid.pixel_grid.n_columns);
  cudaThreadSynchronize();

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
                     sigma_z,
                     sigma_dl,
                     grid);
      cudaThreadSynchronize();

      progress(ib * n_iterations_in_block + it, true);
    }

    output_rho.copy_from_device();
    Output rho_output(grid.pixel_grid.n_columns,
                      grid.pixel_grid.n_rows,
                      grid.n_planes,
                      output_rho.host_ptr);
    output((ib + 1) * n_iterations_in_block, rho_output);
  }

  progress(n_iteration_blocks * n_iterations_in_block, false);
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

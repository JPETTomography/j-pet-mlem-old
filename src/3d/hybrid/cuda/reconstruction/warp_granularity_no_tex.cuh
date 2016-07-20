#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "../reconstruction.h"

#include "common/cuda/kernels.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {
namespace Reconstruction {

__global__ static void reconstruction(const LineSegment* lor_line_segments,
                                      const Pixel* pixels,
                                      const F* pixel_weights,
                                      const Event* events,
                                      const int n_events,
                                      F* output_rho,
                                      F* rho,
                                      F* sensitivity,
                                      const F sigma_z,
                                      const F sigma_dl,
                                      const Grid grid,
                                      const F barrel_length) {

  const Kernel kernel(sigma_z, sigma_dl);
  const auto tid = (blockIdx.x * blockDim.x) + threadIdx.x;
  const auto n_warps_per_block = blockDim.x / WARP_SIZE;
  const auto n_warps = gridDim.x * n_warps_per_block;
  const auto n_chunks = (n_events + n_warps - 1) / n_warps;
  const auto warp_index = tid / WARP_SIZE;

  // --- event loop ----------------------------------------------------------
  for (int chunk_index = 0; chunk_index < n_chunks; ++chunk_index) {
    int event_index =
        chunk_index * n_warps + blockIdx.x * n_warps_per_block + warp_index;
    // check if we are still on the list
    if (event_index >= n_events)
      break;

    const auto event = events[event_index];
    const auto lor = event.lor;
    const auto lor_index = lor.index();
    const auto segment = lor_line_segments[lor_index];
    const auto R = segment.length / 2;
    F denominator = 0;
    const auto n_pixels = event.pixel_info_end - event.pixel_info_begin;
    const auto n_pixel_chunks = (n_pixels + WARP_SIZE + 1) / WARP_SIZE;

    // -- voxel loop - denominator -------------------------------------------
    for (auto pixel_chunk = 0; pixel_chunk < n_pixel_chunks; ++pixel_chunk) {
      const auto pixel_index =
          WARP_SIZE * pixel_chunk + (threadIdx.x & (WARP_SIZE - 1));
      // check if we are still on the list
      if (pixel_index >= n_pixels)
        break;
      const auto info_index = pixel_index + event.pixel_info_begin;
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      const auto center = grid.pixel_grid.center_at(pixel);
      const auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        int voxel_index = grid.index(voxel);
        auto z = grid.center_z_at(voxel);
        auto kernel2d = kernel.normalized(Point2D(event.right, event.up),
                                          event.tan,
                                          event.sec,
                                          R,
                                          barrel_length,
                                          Point2D(z, up));
        if (kernel2d <= 0)  // can happen when pixel is at radius boundary
          continue;
        auto kernel_t = pixel_weight;
        auto weight = kernel2d * kernel_t *  // hybrid of 2D x-y & y-z
                      rho[voxel_index];
        // end of kernel calculation
        denominator +=
            weight *
            sensitivity[grid.pixel_grid.index(Pixel(voxel.x, voxel.y))];
      }
    }  // voxel loop - denominator

    // reduce denominator so all threads now share same value
    Common::GPU::reduce(denominator);

    if (denominator <= 0)  // boundary condition: event at the edge of FOV
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto pixel_chunk = 0; pixel_chunk < n_pixel_chunks; ++pixel_chunk) {
      const auto pixel_index =
          WARP_SIZE * pixel_chunk + (threadIdx.x & (WARP_SIZE - 1));
      // check if we are still on the list
      if (pixel_index >= n_pixels)
        break;
      const auto info_index = pixel_index + event.pixel_info_begin;
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      const auto center = grid.pixel_grid.center_at(pixel);
      const auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        int voxel_index = grid.index(voxel);
        auto z = grid.center_z_at(voxel);
        auto kernel2d = kernel.normalized(Point2D(event.right, event.up),
                                          event.tan,
                                          event.sec,
                                          R,
                                          barrel_length,
                                          Point2D(z, up));
        if (kernel2d <= 0)  // can happen when pixel is at radius boundary
          continue;
        auto kernel_t = pixel_weight;
        auto weight = kernel2d * kernel_t *  // hybrid of 2D x-y & y-z
                      rho[voxel_index];
        // end of kernel calculation

        atomicAdd(&output_rho[voxel_index], weight * inv_denominator);
      }
    }  // voxel loop
  }    // event loop
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

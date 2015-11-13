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
    const auto n_pixels =
        event.last_pixel_info_index - event.first_pixel_info_index;
    const auto n_planes = event.last_plane - event.first_plane;
    const auto n_voxels = n_pixels * n_planes;
    const auto n_voxel_chunks = (n_voxels + WARP_SIZE + 1) / WARP_SIZE;

    // -- voxel loop - denominator -------------------------------------------
    for (auto voxel_chunk = 0; voxel_chunk < n_voxel_chunks; ++voxel_chunk) {
      const auto voxel_index =
          WARP_SIZE * voxel_chunk + (threadIdx.x & (WARP_SIZE - 1));
      // check if we are still on the list
      if (voxel_index >= n_voxels)
        break;
      const auto pixel_index = voxel_index / n_planes;
      const auto plane_index = voxel_index % n_planes;
      const auto& pixel_info =
          pixel_infos[pixel_index + event.first_pixel_info_index];

      auto pixel = pixel_info.pixel;
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);
      auto iz = event.first_plane + plane_index;

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
    }  // voxel loop - denominator

    // reduce denominator so all threads now share same value
    Common::GPU::reduce(denominator);

    if (denominator == 0)
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto voxel_chunk = 0; voxel_chunk < n_voxel_chunks; ++voxel_chunk) {
      const auto voxel_index =
          WARP_SIZE * voxel_chunk + (threadIdx.x & (WARP_SIZE - 1));
      // check if we are still on the list
      if (voxel_index >= n_voxels)
        break;
      const auto pixel_index = voxel_index / n_planes;
      const auto plane_index = voxel_index % n_planes;
      const auto& pixel_info =
          pixel_infos[pixel_index + event.first_pixel_info_index];

      auto pixel = pixel_info.pixel;
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);
      auto iz = event.first_plane + plane_index;

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

      atomicAdd(&output_rho[grid.index(voxel)], weight * inv_denominator);
    }  // voxel loop
  }    // event loop
}

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

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

#if !USE_MULTI_PLANE_WEIGHTS && !USE_3D_SENSITIVITY
texture<F, 3, cudaReadModeElementType> tex_rho;
texture<F, 2, cudaReadModeElementType> tex_sensitivity;
texture<F, 3, cudaReadModeElementType> tex_3d_sensitivity;
#endif

#if USE_MULTI_PLANE_WEIGHTS
namespace Multiplane {
#elif USE_3D_SENSITIVITY
namespace Sensitivity3D {
#endif

__global__ static void reconstruction(const LineSegment* lor_line_segments,
                                      const Pixel* pixels,
                                      const F* pixel_weights,
#if USE_MULTI_PLANE_WEIGHTS
                                      const int n_pixel_infos,
                                      const int n_planes_half,
#elif USE_3D_SENSITIVITY
                                      const int sensitivity_depth,
#endif
                                      const Event* events,
                                      const int n_events,
                                      F* output_rho,
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
#if !USE_MULTI_PLANE_WEIGHTS
      const auto pixel_weight = pixel_weights[info_index];
#endif
      const auto center = grid.pixel_grid.center_at(pixel);
      const auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
#if USE_MULTI_PLANE_WEIGHTS
        const auto abs_plane = compat::abs(iz - n_planes_half);
        const auto pixel_weight =
            pixel_weights[abs_plane * n_pixel_infos + info_index];
#elif USE_3D_SENSITIVITY
        const auto abs_plane = compat::abs(iz - sensitivity_depth);
#endif
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        auto rho = tex3D(tex_rho, voxel.x, voxel.y, voxel.z);
        if (rho == 0)
          continue;
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
                      rho;
        // end of kernel calculation
        denominator += weight *
#if USE_MULTI_PLANE_WEIGHTS || USE_3D_SENSITIVITY
                       tex3D(tex_3d_sensitivity, voxel.x, voxel.y, abs_plane)
#else
                       tex2D(tex_sensitivity, voxel.x, voxel.y)
#endif
            ;
      }
    }  // voxel loop - denominator

    // reduce denominator so all threads now share same value
    Common::GPU::reduce(denominator);

    if (denominator == 0)
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
#if !USE_MULTI_PLANE_WEIGHTS
      const auto pixel_weight = pixel_weights[info_index];
#endif
      const auto center = grid.pixel_grid.center_at(pixel);
      const auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
#if USE_MULTI_PLANE_WEIGHTS
        const auto pixel_weight =
            pixel_weights[compat::abs(iz - n_planes_half) * n_pixel_infos +
                          info_index];
#endif
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
        auto rho = tex3D(tex_rho, voxel.x, voxel.y, voxel.z);
        if (rho == 0)
          continue;
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
                      rho;
        // end of kernel calculation

        int voxel_index = grid.index(voxel);
        atomicAdd(&output_rho[voxel_index], weight * inv_denominator);
      }
    }  // voxel loop
  }    // event loop
}

#if USE_MULTI_PLANE_WEIGHTS || USE_3D_SENSITIVITY
}  // Multiplane || Sensitivity3D
#endif
}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D

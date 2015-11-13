#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

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

}  // Reconstruction
}  // GPU
}  // Hybrid
}  // PET3D
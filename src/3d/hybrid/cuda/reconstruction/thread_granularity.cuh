#include <cuda_runtime.h>

#include "util/cuda/debug.h"  // catches all CUDA errors
#include "util/cuda/memory.h"
#include "util/delegate.h"

#include "../reconstruction.h"

namespace PET3D {
namespace Hybrid {
namespace GPU {
namespace Reconstruction {

texture<F, 3, cudaReadModeElementType> tex_rho;
texture<F, 2, cudaReadModeElementType> tex_sensitivity;

__global__ static void reconstruction(const LineSegment* lor_line_segments,
                                      const Pixel* pixels,
                                      const F* pixel_weights,
                                      const Event* events,
                                      const int n_events,
                                      F* output_rho,
                                      const F sigma_z,
                                      const F sigma_dl,
                                      const Grid grid,
                                      const F barrel_length) {

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
    for (auto info_index = event.pixel_info_begin;
         info_index < event.pixel_info_end;
         ++info_index) {
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
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
                      tex3D(tex_rho, voxel.x, voxel.y, voxel.z);
        // end of kernel calculation
        denominator += weight * tex2D(tex_sensitivity, voxel.x, voxel.y);
      }
    }  // voxel loop - denominator

    if (denominator == 0)
      continue;

    const auto inv_denominator = 1 / denominator;

    // -- voxel loop ---------------------------------------------------------
    for (auto info_index = event.pixel_info_begin;
         info_index < event.pixel_info_end;
         ++info_index) {
      const auto pixel = pixels[info_index];
      const auto pixel_weight = pixel_weights[info_index];
      auto center = grid.pixel_grid.center_at(pixel);
      auto up = segment.projection_relative_middle(center);

      for (int iz = event.plane_begin; iz < event.plane_end; ++iz) {
        // kernel calculation:
        Voxel voxel(pixel.x, pixel.y, iz);
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

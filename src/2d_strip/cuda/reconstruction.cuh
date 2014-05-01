#pragma once

#include <cuda_runtime.h>
#include "../config.h"
#include "../event.h"
#include "reconstruction_methods.cuh"

#define KERNEL_DEBUG 0

__device__ float lookup_table[200][200];

__shared__ float shmem_z_u[128];
__shared__ float shmem_z_d[128];
__shared__ float shmem_delta_l[128];

template <typename T>
__global__ void reconstruction_2d_strip_cuda(gpu_config::GPU_parameters cfg,
                                             event<T>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho) {

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  // future constant memory allocation on host side !FIXME

  __shared__ float inv_c[3];

  if (tid == 0) {

    inv_c[0] = cfg.inv_pow_sigma_z;
    inv_c[1] = cfg.inv_pow_sigma_z;
    inv_c[2] = cfg.inv_pow_sigma_dl;

    for (int i = 0; i < 1; ++i) {

      //      float z_u = event_list[tid + cfg.number_of_blocks *
      //      blockDim.x].z_u;
      //      float z_d = event_list[tid + cfg.number_of_blocks *
      //      blockDim.x].z_d;
      //      float delta_l = event_list[tid + cfg.number_of_blocks *
      //      blockDim.x].dl;

      float z_u = event_list[i].z_u;
      float z_d = event_list[i].z_d;
      float delta_l = event_list[i].dl;
      float acc = 0.f;

#if KERNEL_DEBUG > 0
      printf("INPUT: %d: Z_U:= %f Z_D:=%f D_L:= %f\n", i, z_u, z_d, delta_l);

#endif

      float sqrt_det_correlation_matrix = sqrt(
          cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

      // angle space transformation
      float tn = event_tan(z_u, z_d, cfg.R_distance);
      float y = event_y(delta_l, tn);
      float z = event_z(z_u, z_d, y, tn);

      float angle = atan(tn);

#if KERNEL_DEBUG > 0
      printf("%d: %f %f %f\n", i, tn, y, z);
#endif

      float cos_ = __cosf(angle);

      float sec_ = float(1.0f) / cos_;
      float sec_sq_ = sec_ * sec_;

      float A = (((T(4.0f) / (cos_ * cos_)) * cfg.inv_pow_sigma_dl) +
                 (T(2.0f) * tn * tn * cfg.inv_pow_sigma_z));
      float B = -T(4.0f) * tn * cfg.inv_pow_sigma_z;
      float C = T(2.0f) * cfg.inv_pow_sigma_z;
      float B_2 = (B / T(2.0f)) * (B / T(2.0f));

      float bb_y = bby(A, C, B_2);

      float bb_z = bbz(A, C, B_2);

      int2 center_pixel = pixel_location(y,
                                         z,
                                         cfg.pixel_size,
                                         cfg.pixel_size,
                                         cfg.grid_size_y_,
                                         cfg.grid_size_z_);
#if KERNEL_DEBUG > 0
      printf("A:= %f B:= %f c:= %f B_2:= %f\n", A, B, C, B_2);
      printf("bb_y:= %f bb_z:= %f\n", bb_y, bb_z);
      printf("Center_Pixel y:= %d z:= %d\n", center_pixel.x, center_pixel.y);
#endif

      // bounding box limits for event
      int2 ur =
          make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
      int2 dl =
          make_int2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
      float2 pp;

#if KERNEL_DEBUG > 0
      printf("iz:= %d Limit:= %d \n", dl.y, ur.y);
      printf("iy:= %d Limit:= %d \n", ur.x, dl.x);
#endif

      for (int iz = dl.y; iz < ur.y; ++iz) {
        for (int iy = ur.x; iy < dl.x; ++iy) {

          pp = pixel_center(iy,
                            iz,
                            cfg.pixel_size,
                            cfg.pixel_size,
                            cfg.grid_size_y_,
                            cfg.grid_size_z_);

          if (in_ellipse(A, B, C, y, z, pp)) {

            float sens = sensitivity(
                pp.x, pp.y, cfg.R_distance, cfg.Scentilator_length / 2.0f);

#if KERNEL_DEBUG > 0
            printf("Pixel(%d,%d): %f %f\n", iy, iz, pp.x, pp.y);
#endif

            pp.x -= y;
            pp.y -= z;

#if KERNEL_DEBUG > 0

            printf("Pixel(%d,%d): SUB: %f %f\n", iy, iz, pp.x, pp.y);
#endif
            T event_kernel = calculate_kernel(y,
                                              tn,
                                              sec_,
                                              sec_sq_,
                                              pp,
                                              inv_c,
                                              cfg,
                                              sqrt_det_correlation_matrix) /
                             sens;

#if KERNEL_DEBUG > 0
            printf("KERNEL: %e SENS: %f\n", event_kernel, sens);
#endif
            acc += event_kernel;
          }
        }
      }
#if KERNEL_DEBUG > 0
      printf("ACC: %e\n", acc);
      printf("TN: %f\n",angle);
      printf("COS: %f\n",__cosf(angle));
#endif
    }
  }
}

#pragma once

#include <cuda_runtime.h>
#include "../config.h"
#include "../event.h"
#include "reconstruction_methods.cuh"

#define WARP_GRANULARITY

//#define EVENT_GRANULARITY

// test flag for first event kernel calcualtions
//#define CHECK_KERNEL

#define IMAGE_SPACE_LINEAR_INDEX(Y, Z) (Y * cfg.n_pixels) + Z
#define BUFFOR_LINEAR_INDEX(Y, Z) \
  (blockIdx.x * cfg.n_pixels * cfg.n_pixels) + (Y * cfg.n_pixels) + Z

#define WARP_SIZE 32

#ifdef WARP_GRANULARITY
template <typename T>
__global__ void reconstruction_2d_strip_cuda(gpu_config::GPU_parameters cfg,
                                             soa_event<float>* soa_data,
                                             event<T>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho,
                                             cudaTextureObject_t tex) {

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  float inv_c[3];

  int offset;

  int number_of_blocks = int(
      ceilf(event_list_size / (cfg.number_of_blocks *
                               (cfg.number_of_threads_per_block / WARP_SIZE))));

  int block_size =
      (cfg.number_of_blocks * (cfg.number_of_threads_per_block / WARP_SIZE));

  int thread_warp_index = floorf(threadIdx.x / WARP_SIZE);

  inv_c[0] = cfg.inv_pow_sigma_z;
  inv_c[1] = cfg.inv_pow_sigma_z;
  inv_c[2] = cfg.inv_pow_sigma_dl;

  float sqrt_det_correlation_matrix =
      sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

  for (int i = 0; i < number_of_blocks; ++i) {

    // calculate offset in event memory location for warp
    int warp_id =
        (i * block_size +
         (blockIdx.x * (cfg.number_of_threads_per_block / WARP_SIZE)) +
         thread_warp_index);

    if ((warp_id < event_list_size)) {

      float z_u = soa_data->z_u[warp_id];
      float z_d = soa_data->z_d[warp_id];
      float delta_l = soa_data->dl[warp_id];

      float acc = 0.f;

      float half_grid_size = 0.5f * cfg.grid_size_y_;
      float half_pixel_size = 0.5f * cfg.pixel_size;

      // angle space transformation
      float tn = event_tan(z_u, z_d, cfg.R_distance);
      float y = event_y(delta_l, tn);
      float z = event_z(z_u, z_d, y, tn);

      float angle = atanf(tn);

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

      int2 ur =
          make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
      int2 dl =
          make_int2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
      float2 pp;

      int2 ul =
          make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                    center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));

      int bb_size = ceilf(((ur.y - ul.y) * (dl.x - ur.x)) / WARP_SIZE);

      int2 start_warp = ul;
      int2 tid_pixel;

      for (int k = 0; k < bb_size; ++k) {

        warp_space_pixel(tid_pixel, offset, start_warp, ul, ur, dl, tid);

        pp = pixel_center(tid_pixel.x,
                          tid_pixel.y,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y_,
                          cfg.grid_size_z_,
                          half_grid_size,
                          half_pixel_size);

        if (in_ellipse(A, B, C, y, z, pp)) {

          pp.x -= y;
          pp.y -= z;

          T event_kernel = calculate_kernel<T>(y,
                                               tn,
                                               sec_,
                                               sec_sq_,
                                               pp,
                                               inv_c,
                                               cfg,
                                               sqrt_det_correlation_matrix) /
                           tex2D<float>(tex, tid_pixel.x, tid_pixel.y);

          acc += event_kernel * tex2D<float>(tex, tid_pixel.x, tid_pixel.y) *
                 rho[IMAGE_SPACE_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)];
        }
      }

      for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
        acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
      }

      float inv_acc = 1.0f / acc;

      start_warp = ul;

      for (int k = 0; k < bb_size; ++k) {

        warp_space_pixel(tid_pixel, offset, start_warp, ul, ur, dl, tid);

        pp = pixel_center(tid_pixel.x,
                          tid_pixel.y,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y_,
                          cfg.grid_size_z_,
                          half_grid_size,
                          half_pixel_size);

        if (in_ellipse(A, B, C, y, z, pp)) {

          pp.x -= y;
          pp.y -= z;

          T event_kernel = calculate_kernel<T>(y,
                                               tn,
                                               sec_,
                                               sec_sq_,
                                               pp,
                                               inv_c,
                                               cfg,
                                               sqrt_det_correlation_matrix) /
                           tex2D<float>(tex, tid_pixel.x, tid_pixel.y);

          atomicAdd(
              &image_buffor[BUFFOR_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)],
              (event_kernel *
               rho[IMAGE_SPACE_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)]) *
                  inv_acc * sec_sq_);
        }
      }
    }
  }
}
#endif

#ifdef EVENT_GRANULARITY

template <typename T>
__global__ void reconstruction_2d_strip_cuda(gpu_config::GPU_parameters cfg,
                                             soa_event<float>* soa_data,
                                             event<T>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho,
                                             cudaTextureObject_t tex) {

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  __shared__ float inv_c[3];

  int block_ = int(ceilf(event_list_size / (cfg.number_of_blocks *
                                            cfg.number_of_threads_per_block)));

  if (threadIdx.x == 0) {

    inv_c[0] = cfg.inv_pow_sigma_z;
    inv_c[1] = cfg.inv_pow_sigma_z;
    inv_c[2] = cfg.inv_pow_sigma_dl;
  }

  __syncthreads();

  T sqrt_det_correlation_matrix =
      sqrt(cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_z * cfg.inv_pow_sigma_dl);

  for (int i = 0; i < block_chunk; ++i) {

    if ((i * cfg.number_of_blocks * cfg.number_of_threads_per_block) + tid <
        event_list_size) {

      for (int j = 0; j < 1; ++j) {

#ifdef AOS_ACCESS
        T z_u = event_list[(i * cfg.number_of_blocks *
                            cfg.number_of_threads_per_block) +
                           tid].z_u;
        T z_d = event_list[(i * cfg.number_of_blocks *
                            cfg.number_of_threads_per_block) +
                           tid].z_d;
        T delta_l = event_list[(i * cfg.number_of_blocks *
                                cfg.number_of_threads_per_block) +
                               tid].dl;
#else
        T z_u = soa_data->z_u[(i * cfg.number_of_blocks *
                               cfg.number_of_threads_per_block) +
                              tid];
        T z_d = soa_data->z_d[(i * cfg.number_of_blocks *
                               cfg.number_of_threads_per_block) +
                              tid];
        T delta_l = soa_data->dl[(i * cfg.number_of_blocks *
                                  cfg.number_of_threads_per_block) +
                                 tid];
#endif

#ifdef CHECK_KERNEL
        if (tid == 0 && i == 0) {

          printf("THREAD\n");
          printf("Z_U: %f Z_D: %f DL: %f\n", z_u, z_d, delta_l);
        }
#endif

        T acc = 0.f;

        T half_grid_size = 0.5f * cfg.grid_size_y_;
        T half_pixel_size = 0.5f * cfg.pixel_size;

        // angle space transformation
        T tn = event_tan(z_u, z_d, cfg.R_distance);
        T y = event_y(delta_l, tn);
        T z = event_z(z_u, z_d, y, tn);

        T angle = atanf(tn);

        T cos_ = __cosf(angle);

        T sec_ = float(1.0f) / cos_;
        T sec_sq_ = sec_ * sec_;

        T A = (((T(4.0f) / (cos_ * cos_)) * cfg.inv_pow_sigma_dl) +
               (T(2.0f) * tn * tn * cfg.inv_pow_sigma_z));
        T B = -T(4.0f) * tn * cfg.inv_pow_sigma_z;
        T C = T(2.0f) * cfg.inv_pow_sigma_z;
        T B_2 = (B / T(2.0f)) * (B / T(2.0f));

        T bb_y = bby(A, C, B_2);

        T bb_z = bbz(A, C, B_2);

        int2 center_pixel = pixel_location(y,
                                           z,
                                           cfg.pixel_size,
                                           cfg.pixel_size,
                                           cfg.grid_size_y_,
                                           cfg.grid_size_z_);

        // bounding box limits for event
        int2 ur =
            make_int2(center_pixel.x - pixels_in_line(bb_y, cfg.pixel_size),
                      center_pixel.y + pixels_in_line(bb_z, cfg.pixel_size));
        int2 dl =
            make_int2(center_pixel.x + pixels_in_line(bb_y, cfg.pixel_size),
                      center_pixel.y - pixels_in_line(bb_z, cfg.pixel_size));
        float2 pp;

        int iter = 0;

        for (int iz = dl.y; iz < ur.y; ++iz) {
          for (int iy = ur.x; iy < dl.x; ++iy) {

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y_,
                              cfg.grid_size_z_,
                              half_grid_size,
                              half_pixel_size);

            //            if (tid == 0 && i == 0) {
            //              printf("Pixel[%d,%d]\n", iy, iz);
            //            }

            if (in_ellipse(A, B, C, y, z, pp)) {

              ++iter;

              pp.x -= y;
              pp.y -= z;

              T event_kernel =
                  calculate_kernel<T>(y,
                                      tn,
                                      sec_,
                                      sec_sq_,
                                      pp,
                                      inv_c,
                                      cfg,
                                      sqrt_det_correlation_matrix) /
                  tex2D<float>(tex, iy, iz);

              acc += event_kernel * tex2D<float>(tex, iy, iz) *
                     rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)];

#ifdef CHECK_KERNEL
              if (tid == 0 && i == 0) {

                printf("PIXEL:[%d,%d] PP:[%f,%f] KERNEL: %e\n",
                       iy,
                       iz,
                       pp.x,
                       pp.y,
                       event_kernel);
              }
#endif
            }
          }
        }

        float inv_acc = 1.0f / acc;

        if (tid == 0 && i == 0) {
          printf("ACC: %f\n", acc);
          // printf("Z_s: %d Z_e: %d\n", dl.y, ur.y);
          // printf("Y_s: %d Y_e: %d\n", ur.x, dl.x);
        }

        for (int iz = dl.y; iz < ur.y; ++iz) {
          for (int iy = ur.x; iy < dl.x; ++iy) {

            pp = pixel_center(iy,
                              iz,
                              cfg.pixel_size,
                              cfg.pixel_size,
                              cfg.grid_size_y_,
                              cfg.grid_size_z_,
                              half_grid_size,
                              half_pixel_size);

            if (in_ellipse(A, B, C, y, z, pp)) {

              // if(tid == 0 && i == 0){printf("Pixel[%d,%d]\n",iy,iz);}

              pp.x -= y;
              pp.y -= z;

              T event_kernel =
                  calculate_kernel<T>(y,
                                      tn,
                                      sec_,
                                      sec_sq_,
                                      pp,
                                      inv_c,
                                      cfg,
                                      sqrt_det_correlation_matrix) /
                  tex2D<float>(tex, iy, iz);

              atomicAdd(&image_buffor[BUFFOR_LINEAR_INDEX(iy, iz)],
                        (event_kernel * rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
                            inv_acc * sec_sq_);

              //              if (tid == 0 && i == 0) {

              //                printf("PIXEL[%d,%d], LOCATION: %d DATA: %f\n",
              //                       iy,
              //                       iz,
              //                       BUFFOR_LINEAR_INDEX(iy, iz),
              //                       (event_kernel *
              //                       rho[IMAGE_SPACE_LINEAR_INDEX(iy, iz)]) *
              //                           inv_acc * sec_sq_);
              //              }
            }
          }
        }
#ifdef CHECK_KERNEL
        if (tid == 0 && i == 0) {
          printf("ACC: %e\n", acc);
          // printf("Z_s: %d Z_e: %d\n", dl.y, ur.y);
          // printf("Y_s: %d Y_e: %d\n", ur.x, dl.x);
        }
#endif
      }
    }
  }
}

#endif

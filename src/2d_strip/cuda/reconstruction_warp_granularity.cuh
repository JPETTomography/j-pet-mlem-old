#pragma once

#include <cuda_runtime.h>

template <typename T>
__global__ void reconstruction_2d_strip_cuda(CUDA::Config cfg,
                                             soa_event<float>* soa_data,
                                             Event<T>* event_list,
                                             int event_list_size,
                                             float* image_buffor,
                                             float* rho,
                                             cudaTextureObject_t tex) {
  short sh_mem_pixel_buffor[20 * 512];

  int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

  float inv_c[3];

  inv_c[0] = cfg.inv_pow_sigma_z;
  inv_c[1] = cfg.inv_pow_sigma_z;
  inv_c[2] = cfg.inv_pow_sigma_dl;

  float half_grid_size = 0.5f * cfg.grid_size_y;
  float half_pixel_size = 0.5f * cfg.pixel_size;

  int offset;

  int number_of_blocks = int(
      ceilf(event_list_size / (cfg.number_of_blocks *
                               (cfg.number_of_threads_per_block / WARP_SIZE))));

  int block_size =
      (cfg.number_of_blocks * (cfg.number_of_threads_per_block / WARP_SIZE));

  int thread_warp_index = floorf(threadIdx.x / WARP_SIZE);

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

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("A: %f B: %f C: %f B_2: %f\n", A, B, C, B_2);
      }

      if (tid == 0 && i == 0) {
        printf("bb_y: %f bb_z: %f \n", bb_y, bb_z);
      }
#endif

      int2 center_pixel = pixel_location(y,
                                         z,
                                         cfg.pixel_size,
                                         cfg.pixel_size,
                                         cfg.grid_size_y,
                                         cfg.grid_size_z);

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

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("UR:= %d Limit:= %d \n", ur.x, ur.y);
        printf("DL:= %d Limit:= %d \n", dl.x, dl.y);
        printf("iz:= %d Limit:= %d \n", dl.y, ur.y);
        printf("iy:= %d Limit:= %d \n", ur.x, dl.x);
      }
#endif

      int bb_size = ceilf(((ur.y - ul.y) * (dl.x - ur.x)) / WARP_SIZE);

      int loop_index = 0;

      int2 start_warp = ul;
      int2 tid_pixel;

      for (int k = 0; k < bb_size; ++k) {

        warp_space_pixel(tid_pixel, offset, start_warp, ul, ur, dl, tid);

        pp = pixel_center(tid_pixel.x,
                          tid_pixel.y,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y,
                          cfg.grid_size_z,
                          half_grid_size,
                          half_pixel_size);

        if (in_ellipse(A, B, C, y, z, pp)) {

          pp.x -= y;
          pp.y -= z;

          T event_kernel = main_kernel<T>(y,
                                          tn,
                                          sec_,
                                          sec_sq_,
                                          pp,
                                          inv_c,
                                          cfg,
                                          sqrt_det_correlation_matrix) /
                           tex2D<float>(tex, tid_pixel.x, tid_pixel.y);
#if DEBUG
          if (tid < 32 && i == 0) {
                        printf("TID: %d KERNEL: %e PIXEL: %d %d
                        K:%d\n",tid,event_kernel,tid_pixel.x,tid_pixel.y,k);
          }
#endif
          acc += event_kernel * tex2D<float>(tex, tid_pixel.x, tid_pixel.y) *
                 rho[IMAGE_SPACE_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)];

          sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, loop_index, 0)] =
              tid_pixel.x;
          sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, loop_index, 1)] =
              tid_pixel.y;
          loop_index++;
        }
      }

      for (int xor_iter = 16; xor_iter >= 1; xor_iter /= 2) {
        acc += __shfl_xor(acc, xor_iter, WARP_SIZE);
      }

#if DEBUG
      if (tid == 0 && i == 0) {
        printf("KERNEL: %e\n", acc);
      }
#endif

      float inv_acc = 1.0f / acc;

      for (int k = 0; k < loop_index; ++k) {

        tid_pixel.x = sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, k, 0)];
        tid_pixel.y = sh_mem_pixel_buffor[SH_MEM_INDEX(threadIdx.x, k, 1)];

        pp = pixel_center(tid_pixel.x,
                          tid_pixel.y,
                          cfg.pixel_size,
                          cfg.pixel_size,
                          cfg.grid_size_y,
                          cfg.grid_size_z,
                          half_grid_size,
                          half_pixel_size);

        pp.x -= y;
        pp.y -= z;

        T event_kernel = main_kernel<T>(y,
                                        tn,
                                        sec_,
                                        sec_sq_,
                                        pp,
                                        inv_c,
                                        cfg,
                                        sqrt_det_correlation_matrix) /
                         tex2D<float>(tex, tid_pixel.x, tid_pixel.y);

#if DEBUG
        if (tid < 32 && i == 0 && k < loop_index) {
                    printf("TID: %d KERNEL: %e PIXEL: %d %d
                    K:%d\n",tid,event_kernel,tid_pixel.x,tid_pixel.y,k);
        }
#endif

        atomicAdd(&image_buffor[BUFFOR_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)],
                  (event_kernel *
                   rho[IMAGE_SPACE_LINEAR_INDEX(tid_pixel.x, tid_pixel.y)]) *
                      inv_acc * sec_sq_);
      }
    }
  }
}

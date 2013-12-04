#pragma once

#include <cuda_runtime.h>

#include "data_structures.h"
#include "prng.cuh"
#include "geometry_methods.cuh"

__device__ int lor_iterator(int& id1, int& id2) {

  if (id1 < id2) {
    int temp;
    temp = id2;
    id2 = id1;
    id1 = temp;
  }

  return ((id1 * (id1 + 1)) / 2) + id2;
}

__global__ void gpu_phantom_generation(int x,
                                       int y,
                                       int iteration,
                                       unsigned int* gpu_prng_seed,
                                       Matrix_Element* pixel_data,
                                       int threads,
                                       int pixels_in_row,
                                       float radius,
                                       float h_detector,
                                       float w_detector,
                                       float pixel_size) {

  int tid = ((blockIdx.x * blockDim.x) + threadIdx.x);

  unsigned int seed[4];

  seed[0] = gpu_prng_seed[4 * tid];
  seed[1] = gpu_prng_seed[4 * tid + 1];
  seed[2] = gpu_prng_seed[4 * tid + 2];
  seed[3] = gpu_prng_seed[4 * tid + 3];

  __shared__ Detector_Ring ring;

  Hits hit1;
  Hits hit2;

  float fov_radius = radius / M_SQRT2;
  // FIXME: will fail for large number of detectors
  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    create_detector_ring(h_detector, w_detector, radius, ring);
  }

  __syncthreads();

  int detector1;
  int detector2;

  int i_inner;
  int i_outer;

#pragma unroll
  for (int i = 0; i < iteration; ++i) {

    float rx =
        (x + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;
    float ry =
        (y + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;

    float angle = HybridTaus(seed[0], seed[1], seed[2], seed[3]) * M_PI;

    if (rx > ry) {
      continue;
    }

    // innetr and outer secant for circles
    Secant_Points inner_secant = secant(rx, ry, angle, radius);
    Secant_Points outer_secant = secant(rx, ry, angle, radius + h_detector);

    // hits per detector(if hits = 2 we got pair of detector, else generate
    // new random position and angle)

    i_inner = section(secant_angle(inner_secant.x1, inner_secant.y1),
                      NUMBER_OF_DETECTORS);
    i_outer = section(secant_angle(outer_secant.x1, inner_secant.y1),
                      NUMBER_OF_DETECTORS);

    if (!check_for_hits(i_inner,
                        i_outer,
                        rx,
                        ry,
                        angle,
                        NUMBER_OF_DETECTORS,
                        ring,
                        detector1,
                        hit1)) {
      continue;
    }

    i_inner = section(secant_angle(inner_secant.x2, inner_secant.y2),
                      NUMBER_OF_DETECTORS);
    i_outer = section(secant_angle(outer_secant.x2, inner_secant.y2),
                      NUMBER_OF_DETECTORS);

    if (!check_for_hits(i_inner,
                        i_outer,
                        rx,
                        ry,
                        angle,
                        NUMBER_OF_DETECTORS,
                        ring,
                        detector2,
                        hit2)) {
      continue;
    }

    //    float deposition_depth =
    //        -log(HybridTaus(seed[0], seed[1], seed[2], seed[3])) *
    // inv_unit_prob_;

    // if(threadIdx.x == 1){
    // printf("TID:%d %d %d %d %d %f %f %f %f %f %f %f %f\n ", threadIdx.x, x,
    // y,
    //        detector1, detector2, hit1.p[0].x, hit1.p[0].y, hit1.p[1].x,
    //        hit1.p[1].y, hit2.p[0].x, hit2.p[0].y, hit2.p[1].x, hit2.p[1].y);
    // }
    // if(deposition_depth < (sqrt( (hit[1].x - hit[0].x) * (hit[1].x -
    // hit[0].x) + (hit[1].y - hit[0].y) * (hit[1].x - hit[0].x) ))) {

    atomicAdd(&pixel_data[blockIdx.x].lor[lor_iterator(detector1, detector2)],
              1.0f);
    //}
  }

  gpu_prng_seed[4 * tid] = seed[0];
  gpu_prng_seed[4 * tid + 1] = seed[1];
  gpu_prng_seed[4 * tid + 2] = seed[2];
  gpu_prng_seed[4 * tid + 3] = seed[3];
}

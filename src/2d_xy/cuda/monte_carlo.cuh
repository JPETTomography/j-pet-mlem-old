#pragma once

#include <cuda_runtime.h>

#include "geometry.h"
#include "geometry_methods.cuh"
#include "prng.cuh"

__device__ int lor_iterator(int& id1, int& id2) {

  if (id1 < id2) {
    //    int temp;
    //    temp = id2;
    //    id2 = id1;
    //    id1 = temp;
    id1 ^= id2;
    id2 ^= id1;
    id1 ^= id2;
  }

  return ((id1 * (id1 + 1)) / 2) + id2;
}

__global__ void monte_carlo_kernel(int x,
                                   int y,
                                   int iteration,
                                   int n_detectors,
                                   int tof_n_positions,
                                   unsigned int* gpu_prng_seed,
                                   MatrixElement* pixel_data,
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

  __shared__ DetectorRing ring;

  Hits hit1;
  Hits hit2;

  // float fov_radius = radius / M_SQRT2;
  // FIXME: will fail for large number of detectors
  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    create_detector_ring(h_detector, w_detector, radius, ring);
  }

  __syncthreads();

  int detector1, detector2;
  float depth1, depth2, position;
  Point center;
  SecantPoints inner_secant;
  SecantPoints outer_secant;

#pragma unroll
  for (int i = 0; i < iteration; ++i) {

    center.x =
        (x + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;
    center.y =
        (y + HybridTaus(seed[0], seed[1], seed[2], seed[3])) * pixel_size;

    float angle = HybridTaus(seed[0], seed[1], seed[2], seed[3]) * M_PI;

    if (center.x > center.y) {
      continue;
    }

    // inner and outer secant for circles

    secant(inner_secant,
           outer_secant,
           center.x,
           center.y,
           angle,
           radius,
           radius + h_detector);

    // hits per detector(if hits = 2 we got pair of detector, else generate
    // new random position and angle)

    SecantSections i_inner = secant_sections(inner_secant, NUMBER_OF_DETECTORS);
    SecantSections i_outer = secant_sections(outer_secant, NUMBER_OF_DETECTORS);

    int intersection_flag = 1;

    //    unsigned int t0 = clock();

    intersection_flag = check_for_hits(intersection_flag,
                                       i_inner.ss1,
                                       i_outer.ss1,
                                       center.x,
                                       center.y,
                                       angle,
                                       NUMBER_OF_DETECTORS,
                                       ring,
                                       detector1,
                                       hit1,
                                       seed,
                                       depth1);

    intersection_flag = check_for_hits(intersection_flag,
                                       i_inner.ss2,
                                       i_outer.ss2,
                                       center.x,
                                       center.y,
                                       angle,
                                       NUMBER_OF_DETECTORS,
                                       ring,
                                       detector2,
                                       hit2,
                                       seed,
                                       depth2);

    if (intersection_flag) {
      float length1 = nearest_distance(hit1.p[0], hit1.p[1], center) + depth1;
      float length2 = nearest_distance(hit2.p[0], hit2.p[1], center) + depth2;

      if (detector1 > detector2) {
        position = length1 - length2;
      } else {
        position = length2 - length1;
      }
#if NO_TOF > 0
      atomicAdd(&pixel_data[blockIdx.x].hit[lor_iterator(detector1, detector2)],
                1.0f);
#else
      atomicAdd(&pixel_data[(blockIdx.x * tof_n_positions) +
                            quantize_position(position, 0.01f, tof_n_positions)]
                     .hit[lor_iterator(detector1, detector2)],
                1.0f);

#endif
    }
  }

  gpu_prng_seed[4 * tid] = seed[0];
  gpu_prng_seed[4 * tid + 1] = seed[1];
  gpu_prng_seed[4 * tid + 2] = seed[2];
  gpu_prng_seed[4 * tid + 3] = seed[3];
}

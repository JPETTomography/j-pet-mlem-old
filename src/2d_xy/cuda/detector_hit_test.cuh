#pragma once

#include <cuda_runtime.h>

#include "config.h"
#include "geometry.h"
#include "geometry_methods.cuh"

__global__ void detector_hit_test_kernel(float crx,
                                         float cry,
                                         float cangle,
                                         float radius,
                                         float h_detector,
                                         float w_detector) {

  __shared__ DetectorRing test_ring;

  Hits hit1;
  Hits hit2;

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    create_detector_ring(h_detector, w_detector, radius, test_ring);
  }

  __syncthreads();

  int detector1;
  int detector2;

  if (threadIdx.x == 0 && blockIdx.x == 0) {

    float rx = crx;
    float ry = cry;
    float angle = cangle;

    printf("DATA: rx:= %f ry:= %f angle:= %f\n", rx, ry, angle);

    SecantPoints inner_secant = secant(rx, ry, angle, radius);
    SecantPoints outer_secant = secant(rx, ry, angle, radius + h_detector);

    SecantSections i_inner = secant_sections(inner_secant, NUMBER_OF_DETECTORS);
    SecantSections i_outer = secant_sections(outer_secant, NUMBER_OF_DETECTORS);

    check_for_hits(i_inner.ss1,
                   i_outer.ss2,
                   rx,
                   ry,
                   angle,
                   NUMBER_OF_DETECTORS,
                   test_ring,
                   detector1,
                   hit1);

    check_for_hits(i_inner.ss2,
                   i_outer.ss2,
                   rx,
                   ry,
                   angle,
                   NUMBER_OF_DETECTORS,
                   test_ring,
                   detector2,
                   hit2);

    printf("GPU HIT1: %f %f %f %f\n",
           hit1.p[0].x,
           hit1.p[0].y,
           hit1.p[1].x,
           hit1.p[1].y);

    printf("GPU HIT2: %f %f %f %f\n",
           hit2.p[0].x,
           hit2.p[0].y,
           hit2.p[1].x,
           hit2.p[1].y);
  }
}

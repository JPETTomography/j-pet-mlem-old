#pragma once

#include <cuda_runtime.h>

#include "config.h"
#include "data_structures.h"
#include "geometry_methods.cuh"

__global__ void gpu_detector_hit_test(float crx,
                                      float cry,
                                      float cangle,
                                      float radius,
                                      float h_detector,
                                      float w_detector) {

  __shared__ Detector_Ring test_ring;

  Hits hit1;
  Hits hit2;

  if (threadIdx.x < NUMBER_OF_DETECTORS) {


    create_detector_ring(h_detector, w_detector, radius, test_ring);

  }

  __syncthreads();

  int detector1;
  int detector2;

  int i_inner;
  int i_outer;

  if (threadIdx.x == 0 && blockIdx.x == 0) {

    float rx = crx;
    float ry = cry;
    float angle = cangle;

    printf("DATA: rx:= %f ry:= %f angle:= %f\n", rx, ry, angle);

    Secant_Points inner_secant = secant(rx, ry, angle, radius);
    Secant_Points outer_secant = secant(rx, ry, angle, radius + h_detector);

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
                        test_ring,
                        detector1,
                        hit1)) {
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
                        test_ring,
                        detector2,
                        hit2)) {
    }

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

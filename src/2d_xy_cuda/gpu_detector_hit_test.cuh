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

    Detectors detector_base;

    detector_base.points[0].x =
        (w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[0].y = h_detector / 2.0f;
    detector_base.points[1].x =
        (w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[1].y = -h_detector / 2.0f;
    detector_base.points[2].x =
        (-w_detector / 2.0f) + radius + (h_detector / 2.0f);
    detector_base.points[2].y = -h_detector / 2.0f;
    detector_base.points[3].x =
        (-w_detector / 2.0) + radius + (h_detector / 2.0f);
    detector_base.points[3].y = h_detector / 2.0f;

    test_ring.detector_list[threadIdx.x] = detector_base;

    float angle = 2.0f * M_PI * threadIdx.x / NUMBER_OF_DETECTORS;
    float sin_phi = __sinf(angle);
    float cos_phi = __cosf(angle);

    for (int j = 0; j < 4; ++j) {

      float temp_x = test_ring.detector_list[threadIdx.x].points[j].x;
      float temp_y = test_ring.detector_list[threadIdx.x].points[j].y;

      test_ring.detector_list[threadIdx.x].points[j].x =
          temp_x * cos_phi - temp_y * sin_phi;
      test_ring.detector_list[threadIdx.x].points[j].y =
          temp_x * sin_phi + temp_y * cos_phi;
    }
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

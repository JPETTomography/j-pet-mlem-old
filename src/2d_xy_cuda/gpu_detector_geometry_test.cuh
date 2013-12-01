#pragma once

#include "data_structures.h"
#include "geometry_methods.cuh"

__global__ void gpu_detector_geometry_test(float radius,
                                           float h_detector,
                                           float w_detector,
                                           float pixel_size,
                                           Detector_Ring* cpu_output) {

  __shared__ Detector_Ring test_ring;

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

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    cpu_output->detector_list[threadIdx.x] =
        test_ring.detector_list[threadIdx.x];
  }
}

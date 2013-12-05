#pragma once

#include <cuda_runtime.h>

#include "data_structures.h"
#include "geometry_methods.cuh"

__global__ void gpu_detector_geometry_test(float radius,
                                           float h_detector,
                                           float w_detector,
                                           float pixel_size,
                                           Detector_Ring* cpu_output) {

  __shared__ Detector_Ring test_ring;

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    create_detector_ring(h_detector, w_detector, radius, test_ring);
  }

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    cpu_output->detector_list[threadIdx.x] =
        test_ring.detector_list[threadIdx.x];
  }
}

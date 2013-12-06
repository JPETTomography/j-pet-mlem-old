#pragma once

#include <cuda_runtime.h>

#include "config.h"
#include "geometry.h"
#include "geometry_methods.cuh"

__global__ void detector_geometry_test_kernel(float radius,
                                              float h_detector,
                                              float w_detector,
                                              float pixel_size,
                                              DetectorRing* cpu_output) {

  __shared__ DetectorRing test_ring;

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    create_detector_ring(h_detector, w_detector, radius, test_ring);
  }

  if (threadIdx.x < NUMBER_OF_DETECTORS) {

    cpu_output->detector_list[threadIdx.x] =
        test_ring.detector_list[threadIdx.x];
  }
}

#pragma once

#include <cuda_runtime.h>

#include "2d/barrel/simple_geometry.h"
#include "common/types.h"

namespace Common {
namespace GPU {

using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
using PixelInfo = SimpleGeometry::PixelInfo;

// calculates sensitivity out of given pixel_infos
__global__ void reduce_to_sensitivity(const PixelInfo* pixel_infos,
                                      const size_t n_pixel_infos,
                                      float* output,
                                      const int width);

// invert given values i -> 1/i
__global__ void invert(float* input_output, const size_t size);

}  // GPU
}  // Common

#pragma once

#include "../../geometry/point.h"
#include "../../geometry/pixel.h"

#include "../lor.h"
#include "../simple_geometry.h"
#include "../lm_reconstruction.h"

#include "common/types.h"
#include "util/delegate.h"

namespace PET2D {
namespace Barrel {
/// CUDA optimized subimplementation
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = Barrel::LOR<S>;
using Event = Barrel::LMReconstruction<F, S>::Event;
using SimpleGeometry = Barrel::SimpleGeometry<F, S, Hit>;
using PixelInfo = SimpleGeometry::PixelInfo;

/// \endcond

namespace LMReconstruction {

/// CUDA optimized reconstruction implementation
void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         float sigma,
         int width,
         int height,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, float* rho)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> info);

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D

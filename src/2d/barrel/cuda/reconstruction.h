#pragma once

#include "../../geometry/point.h"
#include "../../geometry/pixel.h"

#include "../lor.h"
#include "../simple_geometry.h"
#include "../reconstruction.h"

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
using Mean = Barrel::Reconstruction<F, S, Hit>::Mean;
using SimpleGeometry = Barrel::SimpleGeometry<F, S, Hit>;
using PixelInfo = SimpleGeometry::PixelInfo;

/// \endcond

namespace Reconstruction {

/// CUDA optimized reconstruction implementation
void run(const SimpleGeometry& geometry,
         const Mean* means,
         int n_means,
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

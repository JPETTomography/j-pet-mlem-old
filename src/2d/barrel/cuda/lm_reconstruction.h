#pragma once

#include "../../geometry/point.h"
#include "../../geometry/pixel.h"
#include "../../geometry/pixel_map.h"

#include "../lor.h"
#include "../geometry_soa.h"
#include "../lm_reconstruction.h"

#include "common/types.h"
#include "util/delegate.h"

namespace PET2D {
namespace Barrel {
/// CUDA optimized subimplementation
namespace GPU {
namespace LMReconstruction {

/// \cond PRIVATE
using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = Barrel::LOR<S>;
using Event = Barrel::LMReconstruction<F, S>::Event;
using Geometry = Barrel::GeometrySOA<F, S>;
using Output = PixelMap<Pixel, F>;
/// \endcond

/// CUDA optimized reconstruction implementation
void run(const Geometry& geometry,
         const Event* events,
         int n_events,
         F sigma,
         int width,
         int height,
         int n_iteration_blocks,
         int n_iterations_in_block,
         util::delegate<void(int iteration, const Output& output)> output,
         util::delegate<void(int completed, bool finished)> progress,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char* device_name)> info);

}  // Reconstruction
}  // GPU
}  // Barrel
}  // PET2D

#pragma once

#include "common/types.h"
#include "util/delegate.h"

#include "2d/barrel/simple_geometry.h"
#include "2d/strip/gausian_kernel.h"
#include "3d/geometry/point.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"

#include "../scanner.h"
#include "../reconstruction.h"

namespace PET3D {
namespace Hybrid {
/// CUDA optimized subimplementation
namespace GPU {
namespace Reconstruction {

/// \cond PRIVATE
using Point = PET3D::Point<F>;
using Voxel = PET3D::Voxel<S>;
using LOR = PET2D::Barrel::LOR<S>;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Kernel = PET2D::Strip::GaussianKernel<F>;
using ReconstructionBase = PET3D::Hybrid::Reconstruction<Scanner, Kernel>;
using Event = ReconstructionBase::FrameEvent;
using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
using PixelInfo = SimpleGeometry::PixelInfo;
using Output = VoxelMap<Voxel, F>;
/// \endcond

/// CUDA optimized reconstruction implementation
void run(const SimpleGeometry& geometry,
         const Event* events,
         int n_events,
         float sigma_z,
         float sigma_dl,
         int width,
         int height,
         int depth,
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
}  // Hybrid
}  // PET3D

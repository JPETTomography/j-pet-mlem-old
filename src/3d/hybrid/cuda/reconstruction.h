#pragma once

#include "common/types.h"
#include "util/delegate.h"

#include "2d/barrel/simple_geometry.h"
#include "2d/strip/gaussian_kernel.h"
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
using Point2D = PET2D::Point<F>;
using Vector2D = PET2D::Vector<F>;
using Voxel = PET3D::Voxel<S>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;
using LineSegment = PET2D::LineSegment<F>;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Kernel = PET2D::Strip::GaussianKernel<F>;
using ReconstructionBase = PET3D::Hybrid::Reconstruction<Scanner, Kernel>;
using Event = ReconstructionBase::FrameEvent;
using SimpleGeometry = PET2D::Barrel::SimpleGeometry<F, S, Hit>;
using PixelInfo = SimpleGeometry::PixelInfo;
using Output = VoxelMap<Voxel, F>;
using Sensitivity = PET2D::PixelMap<Pixel, F>;
using Grid = PET3D::VoxelGrid<F, S>;
/// \endcond

/// CUDA optimized reconstruction implementation
void run(const SimpleGeometry& geometry,
         const Sensitivity& sensitivity,
         const Event* events,
         int n_events,
         float sigma_z,
         float sigma_dl,
         const Grid& grid,
         float barrel_length,
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

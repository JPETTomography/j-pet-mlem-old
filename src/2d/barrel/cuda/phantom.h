#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/detector_ring.h"
#if !__CUDACC__
#include "cmdline.h"
#include "2d/barrel/sparse_matrix.h"
#endif

namespace PET2D {
namespace Barrel {
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using Pixel = PET2D::Pixel<>;
using LOR = Barrel::LOR<>;
using SquareDetector = Barrel::SquareDetector<float>;
using DetectorRing = Barrel::DetectorRing<SquareDetector>;

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D

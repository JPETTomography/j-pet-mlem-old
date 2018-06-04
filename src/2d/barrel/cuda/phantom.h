#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/ring_scanner.h"
#if !__CUDACC__
#include "cmdline.h"
#include "2d/barrel/sparse_matrix.h"
#endif
#include "common/types.h"

namespace PET2D {
namespace Barrel {
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = Barrel::LOR<S>;
using SquareDetector = Barrel::SquareDetector<F>;
using Scanner = Barrel::RingScanner<SquareDetector, S>;

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D

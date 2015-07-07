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

namespace PET2D {
namespace Barrel {
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using Pixel = PET2D::Pixel<int>;
using LOR = Barrel::LOR<int>;
using SquareDetector = Barrel::SquareDetector<float>;
using Scanner = Barrel::RingScanner<SquareDetector, short, 192>;

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D

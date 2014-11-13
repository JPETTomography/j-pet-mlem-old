#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/detector_ring.h"
#include "2d/barrel/model.h"
#if !__CUDACC__
#include "cmdline.h"
#include "2d/barrel/sparse_matrix.h"
#endif

namespace PET2D {
namespace Barrel {
namespace GPU {

#if !__CUDACC__
using OutputMatrix = SparseMatrix<Pixel<>, LOR<>>;
OutputMatrix run_matrix(cmdline::parser& cl);
#endif

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using Pixel = PET2D::Pixel<>;
using LOR = Barrel::LOR<>;
using Event = Barrel::Event<float>;
using SquareDetector = Barrel::SquareDetector<float>;
using DetectorRing = Barrel::DetectorRing<SquareDetector>;
using Model = ScintilatorAccept<float>;

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D

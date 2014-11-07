#pragma once

#include "cmdline.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/sparse_matrix.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
using OutputMatrix = SparseMatrix<Pixel<>, LOR<>>;
OutputMatrix run_matrix(cmdline::parser& cl);
}  // GPU
}  // Barrel
}  // PET2D

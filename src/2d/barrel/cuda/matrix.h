#pragma once

#include "cmdline.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/sparse_matrix.h"

namespace PET2D {
namespace Barrel {
namespace GPU {
using OutputMatrix = SparseMatrix<Pixel<>, LOR<>>;
}  // GPU
}  // Barrel
}  // PET2D

PET2D::Barrel::GPU::OutputMatrix run_gpu_matrix(cmdline::parser& cl);

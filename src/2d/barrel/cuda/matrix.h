#pragma once

#include "cmdline.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/sparse_matrix.h"

typedef SparseMatrix<Pixel<>, LOR<>> OutputMatrix;
OutputMatrix run_gpu_matrix(cmdline::parser& cl);

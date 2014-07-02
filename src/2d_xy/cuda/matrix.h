#pragma once

#include "cmdline.h"
#include "geometry/pixel.h"
#include "2d_xy/lor.h"
#include "2d_xy/sparse_matrix.h"

typedef SparseMatrix<Pixel<>, LOR<>> OutputMatrix;
OutputMatrix run_gpu_matrix(cmdline::parser& cl);

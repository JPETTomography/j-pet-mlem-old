#pragma once

#include "cmdline.h"
#include "geometry/pixel.h"
#include "2d_xy/lor.h"
#include "2d_xy/sparse_matrix.h"

SparseMatrix<Pixel<>, LOR<>> run_gpu(cmdline::parser& cl);

#pragma once

#include "3d/hybrid/scanner.h"
#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "common/model.h"
#include "common/types.h"
#if !__CUDACC__
#include "cmdline.h"
#include "2d/barrel/sparse_matrix.h"
#endif

namespace PET3D {
namespace Hybrid {
/// CUDA optimized subimplementation
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<F>;
using Pixel = PET2D::Pixel<S>;
using LOR = PET2D::Barrel::LOR<S>;
using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Event = Scanner::Event;
using Model = Common::ScintillatorAccept<F>;
#if !__CUDACC__
using OutputMatrix = PET2D::Barrel::SparseMatrix<Pixel, LOR, Hit>;
#endif

/// \endcond

/// CUDA optimized Monte-Carlo implementation
class Matrix {
 public:
  Matrix(const F z,
         const Scanner& scanner,
         int n_threads_per_block,
         int n_blocks,
         F pixel_size,
         F length_scale,
         unsigned int* prng_seed);

  ~Matrix();

  void operator()(const Pixel pixel,  ///< pixel to be processed
                  int n_emissions,    ///< numer of emissions
                  int* pixel_hits     ///<[out] result pixel hits
                  );

#if !__CUDACC__
  static OutputMatrix run(cmdline::parser& cl);
#endif

 private:
  const F z;
  Scanner* gpu_scanner;
  const int n_threads_per_block;
  const int n_blocks;
  const F pixel_size;
  const F length_scale;
  unsigned int* gpu_prng_seed;
  const int pixel_hits_count;
  const int pixel_hits_size;
  int* gpu_pixel_hits;
};

}  // GPU
}  // Barrel
}  // PET2D

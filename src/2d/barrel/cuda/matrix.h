#pragma once

#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/barrel/lor.h"
#include "2d/barrel/square_detector.h"
#include "2d/barrel/ring_scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "common/model.h"
#if !__CUDACC__
#include "cmdline.h"
#include "2d/barrel/sparse_matrix.h"
#endif

namespace PET2D {
namespace Barrel {
/// CUDA optimized subimplementation
namespace GPU {

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using Pixel = PET2D::Pixel<short>;
using LOR = Barrel::LOR<short>;
using Event = Barrel::Event<float>;
using SquareDetector = Barrel::SquareDetector<float>;
using Scanner = Barrel::GenericScanner<SquareDetector, MAX_DETECTORS, short>;
using Model = Common::ScintillatorAccept<float>;
#if !__CUDACC__
using OutputMatrix = Barrel::SparseMatrix<Pixel, LOR, typename LOR::S, int>;
#endif

/// \endcond

/// CUDA optimized Monte-Carlo implementation
class Matrix {
 public:
  Matrix(const Scanner& scanner,
         int n_threads_per_block,
         int n_blocks,
         float pixel_size,
         int n_positions,
         float tof_step,
         float length_scale,
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
  Scanner* gpu_scanner;
  const int n_threads_per_block;
  const int n_blocks;
  const float pixel_size;
  const int n_positions;
  const float tof_step;
  const float length_scale;
  unsigned int* gpu_prng_seed;
  const int pixel_hits_count;
  const int pixel_hits_size;
  int* gpu_pixel_hits;
};

}  // GPU
}  // Barrel
}  // PET2D

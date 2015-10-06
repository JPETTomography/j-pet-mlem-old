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
#include "util/delegate.h"

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

/// CUDA sink for unimplemented versions
template <class ScannerClass> class Matrix {
 public:
  static void run(ScannerClass& scanner,
                  int n_blocks,
                  int n_threads_per_block,
                  int n_emissions,
                  double tof_step,
                  int n_pixels,
                  double s_pixel,
                  double length_scale,
                  util::delegate<void(int, bool)> progress,
                  util::delegate<void(LOR, S, Pixel, Hit)> entry) {
    // unused
    (void)scanner, (void)n_blocks, (void)n_threads_per_block, (void)n_emissions,
        (void)tof_step, (void)n_pixels, (void)s_pixel, (void)length_scale,
        (void)progress, (void)entry;
    throw("GPU does not support this scanner type");
  }
};

/// CUDA optimized Monte-Carlo implementation
template <> class Matrix<Scanner> {
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

  static void run(Scanner& scanner,
                  int n_blocks,
                  int n_threads_per_block,
                  int n_emissions,
                  double z_position,
                  int n_pixels,
                  double s_pixel,
                  double length_scale,
                  util::delegate<void(int, bool)> progress,
                  util::delegate<void(LOR, Pixel, Hit)> entry);

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

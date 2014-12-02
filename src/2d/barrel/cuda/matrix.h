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

/// \cond PRIVATE

using Point = PET2D::Point<float>;
using Pixel = PET2D::Pixel<>;
using LOR = Barrel::LOR<>;
using Event = Barrel::Event<float>;
using SquareDetector = Barrel::SquareDetector<float>;
using DetectorRing = Barrel::DetectorSet<SquareDetector>;
using Model = ScintilatorAccept<float>;
#if !__CUDACC__
using OutputMatrix = Barrel::SparseMatrix<Pixel, LOR>;
#endif

class Matrix {
 public:
  Matrix(const DetectorRing& detector_ring,
         int n_threads_per_block,
         int n_blocks,
         float pixel_size,
         int n_positions,
         float tof_step,
         float length_scale,
         unsigned int* prng_seed);

  ~Matrix();

  void operator()(const Pixel pixel,  //< pixel to be processed
                  int n_emissions,    //< numer of emissions
                  int* pixel_hits     //<[out] result pixel hits
                  );

#if !__CUDACC__
  static OutputMatrix run(cmdline::parser& cl);
#endif

 private:
  DetectorRing* gpu_detector_ring;
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

/// \endcond

}  // GPU
}  // Barrel
}  // PET2D

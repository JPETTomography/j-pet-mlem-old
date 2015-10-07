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
#include "util/delegate.h"
#include "util/random.h"

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

/// \endcond

namespace Matrix {

/// CUDA sink for unimplemented versions
template <class ScannerClass>
void run(ScannerClass& scanner,
         util::random::tausworthe& rng,
         int n_emissions,
         double z_position,
         int n_pixels,
         double s_pixel,
         double length_scale,
         util::delegate<void(int, bool)> progress,
         util::delegate<void(LOR, Pixel, Hit)> entry,
         int device,
         int n_blocks,
         int n_threads_per_block,
         util::delegate<void(const char*)> device_name) {
  // unused
  (void)scanner, (void)rng, (void)n_emissions, (void)z_position, (void)n_pixels,
      (void)s_pixel, (void)length_scale, (void)progress, (void)entry,
      (void)device, (void)n_blocks, (void)n_threads_per_block,
      (void)device_name;
  throw("GPU does not support this scanner type");
}

/// CUDA optimized Monte-Carlo implementation for square detector
template <>
void run<Scanner>(Scanner& scanner,
                  util::random::tausworthe& rng,
                  int n_emissions,
                  double z_position,
                  int n_pixels,
                  double s_pixel,
                  double length_scale,
                  util::delegate<void(int, bool)> progress,
                  util::delegate<void(LOR, Pixel, Hit)> entry,
                  int device,
                  int n_blocks,
                  int n_threads_per_block,
                  util::delegate<void(const char*)> device_name);

}  // Matrix
}  // GPU
}  // Hybrid
}  // PET3D

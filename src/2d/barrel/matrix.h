#pragma once

#include "triangular_pixel_map.h"
#include "sparse_matrix.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// 2D barrel PET system matrix
template <typename PixelType,
          typename LORType,
          typename SType = int,
          typename HitType = int>
class Matrix : public TriangularPixelMap<PixelType, SType, HitType> {
  using Base = TriangularPixelMap<PixelType, SType, HitType>;

 public:
  using Pixel = PixelType;
  using LOR = LORType;
  using S = SType;
  using SS = typename std::make_signed<S>::type;
  using Hit = HitType;
  using SparseMatrix = Barrel::SparseMatrix<PixelType, LORType, SType, HitType>;

  Matrix(S n_pixels_in_row,     ///< number of pixels in each directions
         S n_detectors,         ///< number of detectors stored in the matrix
         S n_tof_positions = 1  ///< number of TOF positions (1 for no TOF)
         )
      : Base(n_pixels_in_row),
        n_detectors_(n_detectors),
        n_tof_positions_(n_tof_positions),
        end_lor_(LOR::end_for_detectors(n_detectors)),
        n_emissions_(0) {
    if (n_detectors % 4)
      throw("number of detectors must be multiple of 4");
    if (n_tof_positions < 1)
      throw("number of TOF positions must be equal or greater 1");
  }

  static LOR begin_lor() { return LOR(); }
  const LOR& end_lor() { return end_lor_; }

  void hit_lor(const LOR& lor, S position, S i_pixel, S hits) {
    (void)lor, (void)position, (void)i_pixel, (void)hits;  // unused
    throw(__PRETTY_FUNCTION__);
  }

  S n_detectors() { return n_detectors_; }
  S n_tof_positions() { return n_tof_positions_; }

  S non_zero_lors() { throw(__PRETTY_FUNCTION__); }

  void add_emissions(Hit n_emissions) { n_emissions_ += n_emissions; }

  S n_emissions() { return n_emissions_; }

  void compact_pixel_index(S i_pixel) { (void)i_pixel; }

  Pixel pixel_at_index(S i_pixel) {
    (void)i_pixel;  // unused
    throw(__PRETTY_FUNCTION__);
  }

  SparseMatrix to_sparse() const { throw(__PRETTY_FUNCTION__); }

 private:
  S n_detectors_;
  S n_tof_positions_;
  LOR end_lor_;
  Hit n_emissions_;
};
}  // Barrel
}  // PET2D

#pragma once

#include "triangular_pixel_map.h"
#include "sparse_matrix.h"

template <typename PixelType,
          typename LORType,
          typename SType = int,
          typename HitType = int>
class Matrix : public TriangularPixelMap<PixelType, SType, HitType> {
  typedef TriangularPixelMap<PixelType, SType, HitType> Super;

 public:
  typedef PixelType Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef ::SparseMatrix<PixelType, LORType, SType, HitType> SparseMatrix;

  /// @param n_pixels_in_row number of pixels in each directions
  /// @param n_detectors     number of detectors stored in the matrix
  Matrix(S n_pixels_in_row, S n_detectors, S n_tof_positions = 1)
      : Super(n_pixels_in_row),
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

  void hit_lor(const LOR& lor __attribute__((unused)),
               S position __attribute__((unused)),
               S i_pixel __attribute__((unused)),
               S hits __attribute__((unused)) = 1) {
    throw(__PRETTY_FUNCTION__);
  }

  S n_detectors() { return n_detectors_; }
  S n_tof_positions() { return n_tof_positions_; }

  S non_zero_lors() { throw(__PRETTY_FUNCTION__); }

  void add_emissions(Hit n_emissions) { n_emissions_ += n_emissions; }

  S n_emissions() { return n_emissions_; }

  void compact_pixel_index(S i_pixel __attribute__((unused))) {}

  SparseMatrix to_sparse() const { throw(__PRETTY_FUNCTION__); }

 private:
  S n_detectors_;
  S n_tof_positions_;
  LOR end_lor_;
  Hit n_emissions_;
};

#pragma once

#include "triangular_pixel_map.h"
#include "sparse_matrix.h"
#include "bstream.h"

template <typename LORType, typename SType = int, typename HitType = int>
class Matrix : public TriangularPixelMap<SType, HitType> {
  typedef TriangularPixelMap<SType, HitType> Super;

 public:
  typedef typename Super::Pixel Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef SparseMatrix<LORType, SType, HitType> Sparse;

  /// @param n_pixels    number of pixels in each directions
  /// @param n_detectors number of detectors stored in the matrix
  Matrix(S n_pixels, S n_detectors)
      : Super(n_pixels),
        n_pixels_(n_pixels),
        n_detectors_(n_detectors),
        end_(LOR::end_for_detectors(n_detectors)),
        n_emissions_(0) {
    if (n_pixels % 2)
      throw("number of pixels must be multiple of 2");
    if (n_detectors % 4)
      throw("number of detectors must be multiple of 4");
  }

  static LOR begin() { return LOR(); }
  const LOR end() { return end_; }

  void add_to_t_matrix(const LOR& lor, S i_pixel) {
    throw(__PRETTY_FUNCTION__);
  }

  S n_detectors() { return n_detectors_; }

  S non_zero_lors() { throw(__PRETTY_FUNCTION__); }

  void increase_n_emissions(Hit n_emissions) { n_emissions_ += n_emissions; }

  S n_emissions() { return n_emissions_; }

  void compact_pixel_index(S i_pixel) {}

  Sparse to_sparse() const { throw(__PRETTY_FUNCTION__); }

 private:
  S n_pixels_;
  S n_detectors_;
  LOR end_;
  Hit n_emissions_;
};

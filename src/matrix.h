#pragma once

#include "triangular_pixel_map.h"
#include "sparse_matrix.h"
#include "bstream.h"

template <typename LORType, typename SType = int, typename HitType = int>
class Matrix : public TriangularPixelMap<SType, HitType> {
  typedef TriangularPixelMap<SType, HitType> Super;

 public:
  typedef LORType LOR;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;
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

  Hit operator()(LOR lor, S x, S y) const { throw(__PRETTY_FUNCTION__); }

  S non_zero_lors() { throw(__PRETTY_FUNCTION__); }

  void increase_n_emissions(Hit n_emissions) { n_emissions_ += n_emissions; }

  S n_emissions() { return n_emissions_; }

  static void read_header(ibstream& in,
                          FileInt& in_is_triangular,
                          FileInt& in_n_pixels,
                          FileInt& in_n_emissions,
                          FileInt& in_n_detectors) {
    throw(__PRETTY_FUNCTION__);
  }

  template <class FileWriter> void output_lor_bitmap(FileWriter& fw, LOR& lor) {
    throw(__PRETTY_FUNCTION__);
  }

  bool output_triangular;

  Sparse to_sparse() const { throw(__PRETTY_FUNCTION__); }

 private:
  S n_pixels_;
  S n_detectors_;
  LOR end_;
  Hit n_emissions_;
};

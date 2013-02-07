/// Sparse system matrix binary file format
/// -----------------------------------------------------
/// uint32_t magic       // 'PETp' (triangular) / 'PETP' (full)
/// uint32_t n_pixels_    // half size [PETp] / full size [PETP]
/// uint32_t n_emissions // per pixel
/// uint32_t n_detectors // regardless of magic
/// while (!eof)
///   uint16_t lor_a, lor_b // pair
///   uint32_t pixel_pair_count
///   for(.. count ..)
///     uint16_t pixel_x, pixel_y // half pixels [PETp] / pixels [PETP]
///     uint32_t pixel_hits

#pragma once

#if _OPENMP
#include <omp.h>
#endif

#include "matrix.h"

/// Provides storage for 1/8 PET system matrix
template <typename PixelType,
          typename LORType,
          typename SType = int,
          typename HitType = int>
class MatrixLORMajor : public Matrix<PixelType, LORType, SType, HitType> {
  typedef Matrix<PixelType, LORType, SType, HitType> Super;

 public:
  typedef PixelType Pixel;
  typedef LORType LOR;
  typedef SType S;
  typedef HitType Hit;
  typedef uint8_t BitmapPixel;
  typedef uint32_t FileInt;
  typedef uint16_t FileHalf;
  typedef Hit* Pixels;
  typedef Pixels* TMatrix;

  /// @param n_pixels    number of pixels in each directions
  /// @param n_detectors number of detectors stored in the matrix
  MatrixLORMajor(S n_pixels, S n_detectors)
      : Super(n_pixels, n_detectors),
        n_detectors_(n_detectors),
        n_2_detectors_(2 * n_detectors),
        n_1_detectors_2_(n_detectors + n_detectors / 2),
        n_1_detectors_4_(n_detectors + n_detectors / 4),
        n_lors_(LOR::end_for_detectors(n_detectors).t_index()) {

    // reserve for all lors
    t_matrix_ = new Pixels[n_lors_]();
  }

  ~MatrixLORMajor() {
    for (auto i = 0; i < n_lors_; ++i) {
      if (t_matrix_[i])
        delete[] t_matrix_[i];
    }
    delete[] t_matrix_;
  }

  void hit_lor(LOR& lor, S i_pixel, S hits = 1) {
    auto i_lor = t_lor_index(lor);
    auto pixels = t_matrix_[i_lor];

    // prealocate pixels for specific lor
    // all are nulls upon startup
    if (!pixels) {
#if _OPENMP
// entering critical section, and check again
// because it may have changed in meantime
#pragma omp critical
      if (!(pixels = t_matrix_[i_lor])) {
        pixels = t_matrix_[i_lor] =
                 new Hit[Super::total_n_pixels_in_triangle()]();
      }
#else
      pixels = t_matrix_[i_lor] =
               new Hit[Super::total_n_pixels_in_triangle()]();
#endif
    }
    pixels[i_pixel] += hits;
  }

  Hit operator()(LOR lor, S x, S y) const {
    bool diag;
    S symmetry;
    auto i_pixel = this->pixel_index(x, y, diag, symmetry);
    auto pixels = t_matrix_[lor_index(lor, symmetry)];
    return pixels ? pixels[i_pixel] * (diag ? 2 : 1) : 0;
  }

  S non_zero_lors() {
    S non_zero_lors = 0;
    for (auto i = 0; i < n_lors_; ++i) {
      if (t_matrix_[i])
        ++non_zero_lors;
    }
    return non_zero_lors;
  }

 private:
  // disable copy contructor
  MatrixLORMajor(const MatrixLORMajor& rhs) {}

  static S t_lor_index(LOR& lor) {
    if (lor.first < lor.second) {
      std::swap(lor.first, lor.second);
    }
    return lor.first * (lor.first + 1) / 2 + lor.second;
  }

  /// Computes LOR index based on given symmetry (1 out 8)
  /// @param lor      detector number pair
  /// @param symmetry number (0..7)
  S lor_index(LOR lor, S symmetry) const {
    if (symmetry & 1) {
      lor.first = (n_2_detectors_ - lor.first) % n_detectors_;
      lor.second = (n_2_detectors_ - lor.second) % n_detectors_;
    }
    if (symmetry & 2) {
      lor.first = (n_1_detectors_2_ - lor.first) % n_detectors_;
      lor.second = (n_1_detectors_2_ - lor.second) % n_detectors_;
    }
    if (symmetry & 4) {
      lor.first = (n_1_detectors_4_ - lor.first) % n_detectors_;
      lor.second = (n_1_detectors_4_ - lor.second) % n_detectors_;
    }
    return t_lor_index(lor);
  }

  LOR end_;
  TMatrix t_matrix_;
  S n_detectors_;
  S n_2_detectors_;
  S n_1_detectors_2_;
  S n_1_detectors_4_;
  S n_lors_;
};

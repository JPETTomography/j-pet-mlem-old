#pragma once

#include <cstdint>

namespace PET2D {
namespace Barrel {

/// Triangular pixel map

/// Provides efficient storage for triangular map of pixels.
/// \image html detector_ring2.pdf.png

template <typename PixelType, typename ValueType> class TriangularPixelMap {
 public:
  using Pixel = PixelType;
  using S = typename Pixel::S;
  using I = typename Pixel::I;
  using Value = ValueType;
  using Values = Value*;

  // reserve for pixel stats
  TriangularPixelMap(S n_pixels_in_row)
      : n_pixels_in_row(n_pixels_in_row),
        n_pixels_in_row_half(n_pixels_in_row / 2),
        total_n_pixels_in_triangle(static_cast<I>(n_pixels_in_row) / 2 *
                                   (n_pixels_in_row / 2 + 1) /
                                   2),
        begin_pixel(),
        end_pixel(Pixel::end_for_n_pixels_in_row(n_pixels_in_row)) {
    if (n_pixels_in_row % 2)
      throw("number of pixels must be multiple of 2");
    values = new Value[total_n_pixels_in_triangle]();
  }

  ~TriangularPixelMap() { delete[] values; }

  /// Computes pixel index and determines symmetry number based on pixel
  /// position
  I pixel_index(Pixel p,     ///< pixel coordinate (0..n_pixels, 0..n_pixels)
                bool& diag,  ///<[out] true if abs(x)==abs(y)
                S& symmetry  ///<[out] symmetry number (0..7)
                ) const {
    // shift so 0,0 is now center
    p.x -= n_pixels_in_row_half;
    p.y -= n_pixels_in_row_half;
    // mirror
    symmetry = 0;
    if (p.x < 0) {
      p.x = -p.x - 1;
      symmetry |= 2;
    }
    if (p.y < 0) {
      p.y = -p.y - 1;
      symmetry |= 1;
    }
    // triangulate
    if (p.x > p.y) {
      std::swap(p.x, p.y);
      symmetry |= 4;
    }
    diag = (p.x == p.y);
    return p.index();
  }

  Value operator[](Pixel& p) const {
    bool diag;
    S symmetry;
    auto i_pixel = pixel_index(p, diag, symmetry);
    return values[i_pixel] * (diag ? 2 : 1);
  }

  Value operator[](Pixel p) const {
    bool diag;
    S symmetry;
    auto i_pixel = pixel_index(p, diag, symmetry);
    return values[i_pixel] * (diag ? 2 : 1);
  }

  void hit(S i_pixel, S hits = 1) { values[i_pixel] += hits; }

  const S n_pixels_in_row;
  const S n_pixels_in_row_half;
  const I total_n_pixels_in_triangle;
  const Pixel begin_pixel;
  const Pixel end_pixel;

 private:
  Values values;
};

}  // Barrel
}  // PET2D

#pragma once

#include <cstdint>

namespace PET2D {
namespace Barrel {

template <typename PixelType, typename SType = int, typename HitType = int>
class TriangularPixelMap {
 public:
  typedef PixelType Pixel;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef Hit* Pixels;
  typedef uint8_t BitmapPixel;

  // reserve for pixel stats
  TriangularPixelMap(S n_pixels_in_row)
      : n_pixels_in_row_(n_pixels_in_row),
        n_pixels_in_row_half_(n_pixels_in_row / 2),
        total_n_pixels_in_triangle_(
            n_pixels_in_row / 2 * (n_pixels_in_row / 2 + 1) / 2),
        end_pixel_(Pixel::end_for_n_pixels_in_row(n_pixels_in_row)) {
    if (n_pixels_in_row % 2)
      throw("number of pixels must be multiple of 2");
    t_hits = new Hit[total_n_pixels_in_triangle_]();
  }

  ~TriangularPixelMap() { delete[] t_hits; }

  static Pixel begin_pixel() { return Pixel(); }
  const Pixel& end_pixel() { return end_pixel_; }

  S n_pixels_in_row() const { return n_pixels_in_row_; }
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  S total_n_pixels_in_triangle() const { return total_n_pixels_in_triangle_; }

  /// Computes pixel index and determines symmetry number based on pixel
  /// position
  /// @param x        pixel x coordinate (0..n_pixels)
  /// @param x        pixel y coordinate (0..n_pixels)
  /// @param diag     outputs true if abs(x)==abs(y)
  /// @param symmetry outputs symmetry number (0..7)
  S pixel_index(Pixel p, bool& diag, S& symmetry) const {
    // shift so 0,0 is now center
    p.x -= n_pixels_in_row_half();
    p.y -= n_pixels_in_row_half();
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

  Hit operator[](Pixel& p) const {
    bool diag;
    S symmetry;
    auto i_pixel = pixel_index(p, diag, symmetry);
    return t_hits[i_pixel] * (diag ? 2 : 1);
  }

  Hit operator[](Pixel p) const {
    bool diag;
    S symmetry;
    auto i_pixel = pixel_index(p, diag, symmetry);
    return t_hits[i_pixel] * (diag ? 2 : 1);
  }

  void hit(S i_pixel, S hits = 1) { t_hits[i_pixel] += hits; }

 private:
  Pixels t_hits;
  S n_pixels_in_row_;
  S n_pixels_in_row_half_;
  S total_n_pixels_in_triangle_;
  Pixel end_pixel_;
};
}  // Barrel
}  // PET2D

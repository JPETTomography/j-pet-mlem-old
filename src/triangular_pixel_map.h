#pragma once

template <typename FType = double, typename SType = int, typename HitType = int>
class TriangularPixelMap {
 public:
  typedef FType F;
  typedef SType S;
  typedef typename std::make_signed<S>::type SS;
  typedef HitType Hit;
  typedef Hit* Pixels;
  typedef uint8_t BitmapPixel;
  typedef typename DetectorRing<F>::LOR LOR;

  // reserve for pixel stats
  TriangularPixelMap(S n_pixels_a)
      : n_pixels_in_row_(n_pixels_a),
        n_pixels_in_row_half_(n_pixels_a / 2),
        total_n_pixels_in_triangle_(n_pixels_a / 2 * (n_pixels_a / 2 + 1) / 2) {
    t_hits_ = new Hit[total_n_pixels_in_triangle_]();
  }

  ~TriangularPixelMap() { delete[] t_hits_; }

  S n_pixels_in_row() const { return n_pixels_in_row_; }
  S n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  S total_n_pixels_in_triangle() const { return total_n_pixels_in_triangle_; }

  static constexpr S t_pixel_index(S x, S y) { return y * (y + 1) / 2 + x; }

  /// Computes pixel index and determines symmetry number based on pixel
  /// position
  /// @param x        pixel x coordinate (0..n_pixels)
  /// @param x        pixel y coordinate (0..n_pixels)
  /// @param diag     outputs true if abs(x)==abs(y)
  /// @param symmetry outputs symmetry number (0..7)
  S pixel_index(SS x, SS y, bool& diag, S& symmetry) const {
    // shift so 0,0 is now center
    x -= n_pixels_in_row_half();
    y -= n_pixels_in_row_half();
    // mirror
    symmetry = 0;
    if (x < 0) {
      x = -x - 1;
      symmetry |= 2;
    }
    if (y < 0) {
      y = -y - 1;
      symmetry |= 1;
    }
    // triangulate
    if (x > y) {
      std::swap(x, y);
      symmetry |= 4;
    }
    diag = (x == y);
    return t_pixel_index(x, y);
  }

  Hit hits(S x, S y) const {
    bool diag;
    S symmetry;
    auto i_pixel = pixel_index(x, y, diag, symmetry);
    return t_hits_[i_pixel] * (diag ? 2 : 1);
  }

  void add_hit(S i_pixel) { ++(t_hits_[i_pixel]); }

  void add_hit(S i_pixel, S h) { t_hits_[i_pixel] += h; }

  template <class FileWriter> void output_bitmap(FileWriter& fw) {
    fw.template write_header<BitmapPixel>(n_pixels_in_row(), n_pixels_in_row());
    Hit pixel_max = 0;
    for (auto y = 0; y < n_pixels_in_row(); ++y) {
      for (auto x = 0; x < n_pixels_in_row(); ++x) {
        pixel_max = std::max(pixel_max, hits(x, y));
      }
    }
    auto gain = static_cast<double>(std::numeric_limits<BitmapPixel>::max()) /
                pixel_max;
    for (auto y = 0; y < n_pixels_in_row(); ++y) {
      BitmapPixel row[n_pixels_in_row()];
      for (auto x = 0; x < n_pixels_in_row(); ++x) {
        row[x] = std::numeric_limits<BitmapPixel>::max() - gain * hits(x, y);
      }
      fw.write_row(row);
    }
  }

 private:
  Pixels t_hits_;
  S n_pixels_in_row_;
  S n_pixels_in_row_half_;
  S total_n_pixels_in_triangle_;
};

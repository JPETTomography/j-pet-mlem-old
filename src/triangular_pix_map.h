#pragma once

template <typename F = double, typename HitType = int>
class TriangularPixelMap {
 public:
  typedef HitType hit_type;
  typedef hit_type* pixels_type;
  typedef uint8_t bitmap_pixel_type;
  typedef typename detector_ring<F>::lor_type lor_type;

  // reserve for pixel stats
  TriangularPixelMap(int n_pixels_a)
      : n_pixels_in_row_(n_pixels_a),
        n_pixels_in_row_half_(n_pixels_a / 2),
        total_n_pixels_in_triangle_(n_pixels_a / 2 * (n_pixels_a / 2 + 1) / 2) {
    t_hits_ = new hit_type[total_n_pixels_in_triangle_]();
  }

  ~TriangularPixelMap() { delete[] t_hits_; }

  int n_pixels_in_row() const { return n_pixels_in_row_; }
  int n_pixels_in_row_half() const { return n_pixels_in_row_half_; }
  int total_n_pixels_in_triangle() const { return total_n_pixels_in_triangle_; }

  static constexpr size_t t_pixel_index(size_t x, size_t y) {
    return y * (y + 1) / 2 + x;
  }

  /// Computes pixel index and determines symmetry number based on pixel
  /// position
  /// @param x        pixel x coordinate (0..n_pixels)
  /// @param x        pixel y coordinate (0..n_pixels)
  /// @param diag     outputs true if abs(x)==abs(y)
  /// @param symmetry outputs symmetry number (0..7)
  size_t pixel_index(ssize_t x, ssize_t y, bool& diag, int& symmetry) const {
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

  hit_type hits(size_t x, size_t y) const {
    bool diag;
    int symmetry;
    auto i_pixel = pixel_index(x, y, diag, symmetry);
    return t_hits_[i_pixel] * (diag ? 2 : 1);
  }

  void add_hit(int i_pixel) { ++(t_hits_[i_pixel]); }

  void add_hit(int i_pixel, int h) { t_hits_[i_pixel] += h; }

  template <class FileWriter> void output_bitmap(FileWriter& fw) {
    fw.template write_header<bitmap_pixel_type>(n_pixels_in_row(),
                                                n_pixels_in_row());
    hit_type pixel_max = 0;
    for (auto y = 0; y < n_pixels_in_row(); ++y) {
      for (auto x = 0; x < n_pixels_in_row(); ++x) {
        pixel_max = std::max(pixel_max, hits(x, y));
      }
    }
    auto gain =
        static_cast<double>(std::numeric_limits<bitmap_pixel_type>::max()) /
        pixel_max;
    for (auto y = 0; y < n_pixels_in_row(); ++y) {
      bitmap_pixel_type row[n_pixels_in_row()];
      for (auto x = 0; x < n_pixels_in_row(); ++x) {
        row[x] =
            std::numeric_limits<bitmap_pixel_type>::max() - gain * hits(x, y);
      }
      fw.write_row(row);
    }
  }

 private:
  pixels_type t_hits_;
  int n_pixels_in_row_;
  int n_pixels_in_row_half_;
  int total_n_pixels_in_triangle_;
};

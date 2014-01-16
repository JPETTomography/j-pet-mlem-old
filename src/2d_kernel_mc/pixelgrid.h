#ifndef PIXELGRID_H
#define PIXELGRID_H

#include <vector>

/**
 * @brief PixelGrid represents a a pixel grid defined by the coordinated
 * of the lower left corner, number o columns and rows and pixel size
 */

template <typename F> class PixelGrid {
 public:
  /// Constructs a grid centered at the origin
  PixelGrid(int n_row, int n_col, F pixel_width, F pixel_height)
      : n_row_(n_row),
        n_col_(n_col),
        pixel_width_(pixel_width),
        pixel_height_(pixel_height),
        ll_z_(-n_row_ * pixel_width / 2.0),
        ll_y_(-n_col_ * pixel_height_ / 2.0),
        values_(n_row_ * n_col_) {}

  PixelGrid(int n_row, int n_col, F pixel_width, F pixel_height, F ll_z, F ll_y)
      : n_row_(n_row),
        n_col_(n_col),
        pixel_width_(pixel_width),
        pixel_height_(pixel_height),
        ll_z_(ll_z),
        ll_y_(ll_y),
        values_(n_row_ * n_col_) {}

  int index(int col, int row) const { return row * n_col_ + col; }
  F operator()(int col, int row) const { return values_[index(col, row)]; }
  int col(int z) { return static_cast<int>(floor((z - ll_z_) / pixel_width_)); }
  int row(int y) {
    return static_cast<int>(floor((y - ll_y_) / pixel_height_));
  };

  void add(F z, F y, F weight) { values_[index(col(z), row(y))] += weight; }

 private:
  int n_row_;
  int n_col_;
  F ll_z_;
  F ll_y_;
  F pixel_width_;
  F pixel_height_;

  std::vector<F> values_;
};

#endif  // PIXELGRID_H

#pragma once

#if !__CUDACC__
#include <iostream>
#include "util/bstream.h"
#endif

#include "2d/geometry/point.h"
#include "2d/geometry/vector.h"
#include "2d/geometry/pixel.h"

namespace PET2D {

/// 2D pixel grid description

/// 2D pixel grid description, without actual pixel storage
template <typename FType, typename SType> class PixelGrid {
 public:
  using F = FType;
  using S = SType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Pixel = PET2D::Pixel<S>;

  PixelGrid(S n_columns, S n_rows, F pixel_size, const Point& lower_left)
      : n_columns(n_columns),
        n_rows(n_rows),
        pixel_size(pixel_size),
        lower_left(lower_left),
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)),
        n_pixels(n_columns * n_rows) {}

  const S n_columns;  // n_x
  const S n_rows;     // n_y
  const F pixel_size;
  const Point lower_left;
  const Point lower_left_center;
  const int n_pixels;

  int index(S column, S row) const { return column + n_columns * row; }
  int index(const Pixel& pixel) const { return index(pixel.x, pixel.y); }

  Point lower_left_at(S column, S row) const {
    Vector displacement(column * pixel_size, row * pixel_size);
    return lower_left + displacement;
  }

  Point center_at(S column, S row) const {
    Vector displacement(column * pixel_size, row * pixel_size);
    return lower_left_center + displacement;
  }

  Pixel pixel_at(Point p) const {
    Vector v = p - lower_left;
    S column = static_cast<S>(floor(v.x / pixel_size));
    S row = static_cast<S>(floor(v.y / pixel_size));
    return Pixel(column, row);
  }

  bool contains(Pixel pixel) const {
    return pixel.x >= 0 && pixel.y >= 0 && pixel.x < n_columns &&
           pixel.y < n_rows;
  }

#if !__CUDACC__
  PixelGrid(std::istream& in)
      : n_columns(util::read<S>(in)),
        n_rows(util::read<S>(in)),
        pixel_size(util::read<S>(in)),
        lower_left(in),
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)),
        n_pixels(n_columns * n_rows) {}

  friend util::obstream& operator<<(util::obstream& out, const PixelGrid& pg) {
    out << pg.n_columns << pg.n_rows << pg.pixel_size << pg.lower_left;
    return out;
  }
#endif
};
}

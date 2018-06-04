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
////
/// 2D pixel grid description, without actual pixel storage
template <typename FType, typename SType> class PixelGrid {
 public:
  using F = FType;
  using S = SType;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Pixel = PET2D::Pixel<S>;

  _ PixelGrid(S n_columns, S n_rows, F pixel_size, const Point& lower_left)
      : n_columns(n_columns),
        n_rows(n_rows),
        pixel_size(pixel_size),
        lower_left(lower_left),
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)),
        n_pixels(n_columns * n_rows) {}

  _ PixelGrid(S n_columns, S n_rows, F pixel_size)
      : PixelGrid(
            n_columns,
            n_rows,
            pixel_size,
            Point(-pixel_size * n_columns / 2, -pixel_size * n_rows / 2)) {}

  const S n_columns;  // n_x
  const S n_rows;     // n_y
  const F pixel_size;
  const Point lower_left;
  const Point lower_left_center;
  const int n_pixels;

  _ int index(Pixel pixel) const { return pixel.x + n_columns * pixel.y; }

  _ Point lower_left_at(Pixel pixel) const {
    Vector displacement(pixel.x * pixel_size, pixel.y * pixel_size);
    return lower_left + displacement;
  }

  _ Point center_at(Pixel pixel) const {
    Vector displacement(pixel.x * pixel_size, pixel.y * pixel_size);
    return lower_left_center + displacement;
  }

  _ Pixel pixel_at(Point p) const {
    Vector v = p - lower_left;
    S column = static_cast<S>(compat::floor(v.x / pixel_size));
    S row = static_cast<S>(compat::floor(v.y / pixel_size));
    return Pixel(column, row);
  }

  _ bool contains(Pixel pixel) const {
    return pixel.x >= 0 && pixel.x < n_columns &&  //
           pixel.y >= 0 && pixel.y < n_rows;       //
  }

#if !__CUDACC__
  /// Construct pixel grid from stream
  PixelGrid(std::istream& in)
      : n_columns(util::read<S>(in)),
        n_rows(util::read<S>(in)),
        pixel_size(util::read<F>(in)),
        lower_left(in),
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)),
        n_pixels(n_columns * n_rows) {}

  /// Construct pixel grid from binary stream
  PixelGrid(util::ibstream& in)
      : n_columns(in.read<S>()),
        n_rows(in.read<S>()),
        pixel_size(in.read<F>()),
        lower_left(in),
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)),
        n_pixels(n_columns * n_rows) {}

  /// Output pixel grid to text stream
  friend std::ostream& operator<<(std::ostream& out, const PixelGrid& pg) {
    out << pg.n_columns << ' ' << pg.n_rows << ' ' << pg.pixel_size << ' '
        << pg.lower_left;
    return out;
  }

  /// Output pixel grid to binary stream
  friend util::obstream& operator<<(util::obstream& out, const PixelGrid& pg) {
    out << pg.n_columns << pg.n_rows << pg.pixel_size << pg.lower_left;
    return out;
  }
#endif
};
}

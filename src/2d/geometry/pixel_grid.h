#ifndef PIXEL_GRID
#define PIXEL_GRID

#include <iostream>

#include "2d/geometry/point.h"
#include "2d/geometry/vector.h"
#include "2d/geometry/pixel.h"

namespace PET2D {

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
        lower_left_center(lower_left + Vector(pixel_size / 2, pixel_size / 2)) {
  }

  const S n_columns;  // n_x
  const S n_rows;     // n_y
  const F pixel_size;
  const Point lower_left;
  const Point lower_left_center;

  Point lower_left_at(S column, S row) const {
    Vector displacement(column * pixel_size, row * pixel_size);
    return lower_left + displacement;
  };

  Point center_at(S column, S row) const {
    Vector displacement(column * pixel_size, row * pixel_size);
    return lower_left_center + displacement;
  };

  Pixel pixel_at(Point p) const {
    Vector v = p - lower_left;
    S column = static_cast<S>(floor(v.x / pixel_size));
    S row = static_cast<S>(floor(v.y / pixel_size));
    return Pixel(column, row);
  }

  std::ostream& write(std::ostream& out) {
    out.write((const char*)&n_columns, sizeof(S));
    out.write((const char*)&n_rows, sizeof(S));
    out.write((const char*)&pixel_size, sizeof(F));
    out.write((const char*)&lower_left.x, sizeof(F));
    out.write((const char*)&lower_left.y, sizeof(F));
    return out;
  }

  static PixelGrid read(std::istream& in) {
    S n_columns, n_rows;
    F pixel_size;
    Point lower_left;
    in.read((char*)&n_columns, sizeof(S));
    in.read((char*)&n_rows, sizeof(S));
    in.read((char*)&pixel_size, sizeof(F));
    in.read((char*)&lower_left.x, sizeof(F));
    in.read((char*)&lower_left.y, sizeof(F));

    return PixelGrid(n_columns, n_rows, pixel_size, lower_left);
  }
};
}
#endif  // PIXEL_GRID

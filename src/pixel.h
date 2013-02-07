#pragma once

#include <iostream>
#include <algorithm>

template <typename SType = int> class Pixel {
 public:
  typedef SType S;

  Pixel(S x, S y) : x(x), y(y) {}

  // default constructor
  Pixel() : Pixel(static_cast<S>(0), static_cast<S>(0)) {}

  S x, y;

  constexpr S index() { return y * (y + 1) / 2 + x; }

  Pixel& operator++() {
    if (++x > y) {
      y++;
      x = 0;
    }
    return *this;
  }

  static constexpr Pixel end_for_n_pixels_in_row(S pixels_in_row) {
    return Pixel(0, pixels_in_row);
  }

  bool operator!=(const Pixel& p) const { return x != p.x || y != p.y; }

  bool operator==(const Pixel& p) const { return x == p.x && y == p.y; }

  bool operator<(const Pixel& p) const {
    return y < p.y || (y == p.y && x < p.x);
  }
};

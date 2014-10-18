#pragma once

#include "util/cuda/compat.h"

namespace PET2D {

template <typename SType = int> class Pixel {
 public:
  typedef SType S;

  _ Pixel(S x, S y) : x(x), y(y) {}

  // default constructor
  _ Pixel() : x(0), y(0) {}

  S x, y;

  _ const S index() const { return y * (y + 1) / 2 + x; }

  _ const S index(S width) const { return y * width + x; }

  _ Pixel& operator++() {
    if (++x > y) {
      y++;
      x = 0;
    }
    return *this;
  }

  static const Pixel end_for_n_pixels_in_row(S pixels_in_row) {
    return Pixel(0, pixels_in_row);
  }

  _ bool operator!=(const Pixel& p) const { return x != p.x || y != p.y; }

  _ bool operator==(const Pixel& p) const { return x == p.x && y == p.y; }

  _ bool operator<(const Pixel& p) const {
    return y < p.y || (y == p.y && x < p.x);
  }

  _ void clamp(const Pixel& tl, const Pixel& br) {
    x = compat::min(br.x, compat::max(tl.x, x));
    y = compat::min(br.y, compat::max(tl.y, y));
  }
};
}  // PET2D

#ifdef TEST_CASE
namespace Catch {
template <typename SType> struct StringMaker</**/ ::Pixel<SType>> {
  static std::string convert(const ::Pixel<SType>& p) {
    std::ostringstream oss;
    oss << "(" << p.x << ", " << p.y << ")";
    return oss.str();
  }
};
}
#endif

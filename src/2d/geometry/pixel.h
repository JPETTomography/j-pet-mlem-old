#pragma once

#include "util/cuda/compat.h"
#include "util/read.h"

namespace PET2D {

/// Discreete coordinates pixel
template <typename SType = int> class Pixel {
 public:
  using S = SType;

  _ Pixel(S x, S y) : x(x), y(y) {}

  // default constructor
  _ Pixel() : x(0), y(0) {}

  S x, y;

#if !__CUDACC__
  /// constructs Pixel from stream
  Pixel(std::istream& in) : x(util::read<S>(in)), y(util::read<S>(in)) {}
#endif

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
template <typename SType> struct StringMaker<PET2D::Pixel<SType>> {
  static std::string convert(const PET2D::Pixel<SType>& p) {
    std::ostringstream oss;
    oss << "(" << p.x << ", " << p.y << ")";
    return oss.str();
  }
};
}
#endif

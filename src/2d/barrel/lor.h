#pragma once

#include "util/cuda/compat.h"
#include "util/read.h"

namespace PET2D {
namespace Barrel {

/// Line of Response
template <typename SType> class LOR {
 public:
  using S = SType;

  _ LOR(S first, S second)
      : first(compat::max(first, second)),  // first is always greater
        second(compat::min(first, second)) {}
  _ LOR() = default;

  S first, second;

#if !__CUDACC__
  /// constructs Pixel from stream
  LOR(std::istream& in) : first(util::read<S>(in)), second(util::read<S>(in)) {}
#endif

  _ int index() const { return first * (first + 1) / 2 + second; }

  _ int index(S width) const { return first * width + second; }

  _ LOR& operator++() {
    if (++second > first) {
      first++;
      second = 0;
    }
    return *this;
  }

  static const LOR end_for_detectors(S n_detectors) {
    return LOR(n_detectors, 0);
  }

  _ bool operator==(const LOR& lor) const {
    return second == lor.second && first == lor.first;
  }

  _ bool operator<(const LOR& lor) const {
    return first < lor.first || (first == lor.first && second < lor.second);
  }

  _ bool operator!=(const LOR& lor) const { return !operator==(lor); }

  _ bool operator>(const LOR& lor) const {
    return !(*this < lor) && !(*this == lor);
  }
};
}  // Barrel
}  // PET2D

#ifdef TEST_CASE
namespace Catch {
template <typename SType> struct StringMaker<PET2D::Barrel::LOR<SType>> {
  static std::string convert(const PET2D::Barrel::LOR<SType>& lor) {
    std::ostringstream oss;
    oss << "<" << lor.second << ", " << lor.second << ">";
    return oss.str();
  }
};
}
#endif

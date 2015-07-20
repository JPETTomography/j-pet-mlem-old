#pragma once

#include <type_traits>

#include "util/cuda/compat.h"

#if !__CUDACC__
#include "util/read.h"
#include "util/bstream.h"
#endif

namespace PET2D {
namespace Barrel {

/// Line of Response
template <typename SType> class LOR {
 public:
  using S = SType;
  using Size = typename std::common_type<S, int>::type;

  _ LOR(S first, S second)
      : first(compat::max(first, second)),  // first is always greater
        second(compat::min(first, second)) {}
  _ LOR() = default;

  S first, second;

#if !__CUDACC__
  /// Constructs LOR from stream
  LOR(std::istream& in) : first(util::read<S>(in)), second(util::read<S>(in)) {}

  /// Constructs LOR from binary stream
  LOR(util::ibstream& in) : first(in.read<S>()), second(in.read<S>()) {}
#endif

  _ Size index() const {
    return static_cast<Size>(first) * (first + 1) / 2 + second;
  }

  _ Size index(S width) const {
    return static_cast<Size>(first) * width + second;
  }

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

  _ bool operator!=(const LOR& lor) const { return !operator==(lor); }

  _ bool operator<(const LOR& lor) const {
    return first < lor.first || (first == lor.first && second < lor.second);
  }

  _ bool operator>(const LOR& lor) const {
    return first > lor.first || (first == lor.first && second > lor.second);
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

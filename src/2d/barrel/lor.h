#pragma once

#include <iostream>
#include <algorithm>

namespace PET2D {
namespace Barrel {

/// Line of Response
template <typename SType = int> class LOR : public std::pair<SType, SType> {
 public:
  using S = SType;

  LOR() : std::pair<S, S>(static_cast<S>(0), static_cast<S>(0)) {}

  LOR(S first, S second) : std::pair<S, S>(first, second) {
    if (this->first < this->second) {
      std::swap(this->first, this->second);
    }
  }

  const S index() const {
    return (this->first * (this->first + 1)) / 2 + this->second;
  }

  LOR& operator++() {
    if (++this->second > this->first) {
      this->first++;
      this->second = 0;
    }
    return *this;
  }

  static const LOR end_for_detectors(S n_detectors) {
    return LOR(n_detectors, 0);
  }

  friend std::ostream& operator<<(std::ostream& out, const LOR& lor) {
    return out << '(' << lor.first << ", " << lor.second << ')';
  }
};
}  // Barrel
}  // PET2D

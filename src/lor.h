#pragma once

#include <iostream>
#include <algorithm>

template <typename SType = int> class LOR : public std::pair<SType, SType> {
 public:
  typedef SType S;

  LOR() : std::pair<S, S>(static_cast<S>(0), static_cast<S>(0)) {}

  LOR(S first, S second) : std::pair<S, S>(first, second) {
    if (this->first < this->second) {
      std::swap(this->first, this->second);
    }
  }

  constexpr S index() const {
    return (this->first * (this->first + 1)) / 2 + this->second;
  }

  LOR& operator++() {
    if (++this->second > this->first) {
      this->first++;
      this->second = 0;
    }
    return *this;
  }

  static constexpr LOR end_for_detectors(S n_detectors) {
    return LOR(n_detectors, 0);
  }

  friend std::ostream& operator<<(std::ostream& out __attribute__((unused)),
                                  const LOR& lor __attribute__((unused))) {
    // FIXME: implement me!
    throw(__PRETTY_FUNCTION__);
  }

  /// Computes LOR based on given symmetry (1 out 8)
  /// @param symmetry    number (0..7)
  /// @param n_detectors number of detectors
  LOR symmetric(S symmetry, S n_detectors) const {
    S n_2_detectors = 2 * n_detectors;
    S n_detectors_2 = n_detectors / 2;
    S n_detectors_4 = n_detectors / 4;
    S n_1_detectors_2 = n_detectors + n_detectors_2;
    S n_1_detectors_4 = n_detectors + n_detectors_4;
    LOR lor(*this);
    if (symmetry & 1) {
      lor.first = (n_2_detectors - lor.first) % n_detectors;
      lor.second = (n_2_detectors - lor.second) % n_detectors;
    }
    if (symmetry & 2) {
      lor.first = (n_1_detectors_2 - lor.first) % n_detectors;
      lor.second = (n_1_detectors_2 - lor.second) % n_detectors;
    }
    if (symmetry & 4) {
      lor.first = (n_1_detectors_4 - lor.first) % n_detectors;
      lor.second = (n_1_detectors_4 - lor.second) % n_detectors;
    }
    return lor;
  }
};

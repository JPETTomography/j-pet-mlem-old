#pragma once

#include "cuda/compat.h"

/// Various utility classes and functions

namespace util {

/// Custom random generators and drop-in `std` replacements

namespace random {

/// \a Tausworthe random number generator
class tausworthe {
 public:
  typedef unsigned int result_type;
  typedef long seed_type;
  typedef unsigned int state_type[4];

  /// Minimum value returned by random number generator
  _ static constexpr result_type min() { return 0; }
  /// Maximum value returned by random number generator
  _ static constexpr result_type max() {
    return compat::numeric_max<result_type>();
  }

  /// Creates new random number generator with given seed.
  tausworthe(seed_type seed = 121245) { this->seed(seed); }

  /// Creates new random number generator with seed taken from given memory.
  _ tausworthe(const state_type state) { load(state); }

  /// Creates new random number generator with other RNG.
  tausworthe(tausworthe& rng) {
    for (size_t i = 0; i < sizeof(state) / sizeof(*state); ++i) {
      state[i] = rng();
    }
  }

  /// Loads generator state
  _ void load(const state_type new_state) {
    for (size_t i = 0; i < sizeof(state) / sizeof(*state); ++i) {
      state[i] = new_state[i];
    }
  }

  /// Save generator state
  _ void save(state_type out_state) const {
    for (size_t i = 0; i < sizeof(state) / sizeof(*state); ++i) {
      out_state[i] = state[i];
    }
  }

  /// Compare with other RNG comparing state.
  bool operator==(tausworthe& other) const {
    for (size_t i = 0; i < sizeof(state) / sizeof(*state); ++i) {
      if (other.state[i] != state[i])
        return false;
    }
    return true;
  }

  /// Returns random number
  _ result_type operator()() {
    taus_step(state[0], 13, 19, 12, 4294967294u);
    taus_step(state[1], 2, 25, 4, 4294967288u);
    taus_step(state[2], 3, 11, 17, 4294967280u);
    LCG_step(state[3], 1664525u, 1013904223u);
    return state[0] ^ state[1] ^ state[2] ^ state[3];
  }

  /// Randomizes generator with given seed value
  void seed(seed_type seed) {
#if !_WIN32
    srand48(seed);
#else
    srand(seed);
#endif
    for (size_t i = 0; i < sizeof(state) / sizeof(*state); ++i) {
      result_type r;
#if !_WIN32
      while ((r = static_cast<result_type>(lrand48())) < 128)
#else
      while ((r = static_cast<result_type>(rand())) < 128)
#endif
        ;
      state[i] = r;
    }
  }

 private:
  state_type state;

  template <typename T> _ void taus_step(T& z, int S1, int S2, int S3, T M) {
    unsigned b = (((z << S1) ^ z) >> S2);
    z = (((z & M) << S3) ^ b);
  }

  template <typename T> _ void LCG_step(T& z, T A, T C) { z = (A * z + C); }
};

/// Uniform real distribution for given range and \a RNG
template <typename FType> class uniform_real_distribution {
 public:
  typedef FType result_type;

  /// Creates new distribution with given [a, b) range
  _ uniform_real_distribution(result_type a = 0,  ///< minimum value
                              result_type b = 1   ///< maxumim value
                              )
      : size_(b - a), offset_(a) {}

  /// Returns value from given range using generator
  template <class RNG> _ result_type operator()(RNG& rng) {
    return rng() * size() * scale<RNG>() + offset() -
           static_cast<result_type>(RNG::min()) / range<RNG>();
  }

  /// Return distribution range size
  _ result_type size() const { return size_; }
  /// Return distribution range offset
  _ result_type offset() const { return offset_; }

  template <class RNG> _ static result_type range() {
    return static_cast<result_type>(RNG::max()) -
           static_cast<result_type>(RNG::min());
  }

  template <class RNG> _ static result_type scale() {
    return static_cast<result_type>(1) / range<RNG>();
  }

 private:
  result_type size_, offset_;
};
}  // random
}  // util

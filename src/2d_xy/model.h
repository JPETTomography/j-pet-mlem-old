#pragma once

#include "util/random.h"

/// Represents model which always produces a decay

template <typename FType = double> class AlwaysAccept {
 public:
  typedef FType F;

  AlwaysAccept() {}
  template <class RandomGenerator> bool operator()(RandomGenerator&, F) {
    return true;
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator&) {
    return static_cast<F>(0);
  }

  static F max_bias() { return static_cast<F>(0); }
};

/// Represents model of scintilator where CDF of decay is given by:
/// F = 1-e^(-scale * length), where for length l we get:
/// l = 0:           F(l) = 0
/// l = 1/2 * scale: F(l) = 1 - 1/sqrt(e)
/// l = 1/scale:     F(l) = 1 - 1/e
/// l = 2/scale:     F(l) = 1 - 1/e^2

template <typename FType = double> class ScintilatorAccept {
 public:
  typedef FType F;

  ScintilatorAccept(F scale)
      : one_dis(static_cast<F>(0), static_cast<F>(1)),
        scale(scale),
        inv_scale(static_cast<F>(1) / scale) {}

  template <class RandomGenerator>
  bool operator()(RandomGenerator& gen, F length) {
    return one_dis(gen) >= exp(-length * inv_scale);
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator& gen) {
    auto r = one_dis(gen);
    return -log(r) * scale;
  }

  static F max_bias() { return static_cast<F>(0); }

 private:
  uniform_real_distribution<F> one_dis;
  F scale;
  F inv_scale;
};

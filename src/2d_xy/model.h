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

/// Represents model of scintilator where probability of decay is given by:
/// $p = 1-e^{-scale * length}$, where for:
/// $length = 0$, probability is equal to $0$
/// $length = 1/2*scale$, probability is equal to $1-1/\sqrt{e}$
/// $length = 1/scale$, probability is equal to $1-1/e$
/// $length = 2/scale$, probability is equal to $1-1/e^2$

template <typename FType = double> class ScintilatorAccept {
 public:
  typedef FType F;

  ScintilatorAccept(F scale)
      : one_dis_(static_cast<F>(0), static_cast<F>(1)),
        scale_(scale),
        inv_scale_(static_cast<F>(1) / scale) {}

  template <class RandomGenerator>
  bool operator()(RandomGenerator& gen, F length) {
    return one_dis_(gen) >= exp(-length * inv_scale_);
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator& gen) {
    auto r = one_dis_(gen);
    return -log(r) * scale_;
  }

  static F max_bias() { return static_cast<F>(0); }

 private:
  uniform_real_distribution<F> one_dis_;
  F scale_;
  F inv_scale_;
};

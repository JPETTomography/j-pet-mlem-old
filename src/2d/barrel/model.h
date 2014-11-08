#pragma once

#include "util/random.h"
#include "util/cuda/compat.h"

namespace PET2D {
namespace Barrel {

/// Model which always produces a decay
template <typename FType = double> class AlwaysAccept {
 public:
  using F = FType;

  AlwaysAccept() {}
  template <class RandomGenerator> _ bool operator()(RandomGenerator&, F) {
    return true;
  }

  template <typename RandomGenerator> _ F deposition_depth(RandomGenerator&) {
    return static_cast<F>(0);
  }

  _ static F max_bias() { return static_cast<F>(0); }
};

/// Model of scintilator acceptance

/// Represents model of scintilator where CDF of decay is given by:
/// F = 1-e^(-scale * length), where for length l we get:
/// l = 0:           F(l) = 0
/// l = 1/2 * scale: F(l) = 1 - 1/sqrt(e)
/// l = 1/scale:     F(l) = 1 - 1/e
/// l = 2/scale:     F(l) = 1 - 1/e^2

template <typename FType = double> class ScintilatorAccept {
 public:
  using F = FType;

  _ ScintilatorAccept(F scale)
      : one_dis(static_cast<F>(0), static_cast<F>(1)),
        scale(scale),
        inv_scale(static_cast<F>(1) / scale) {}

  template <class RandomGenerator>
  _ bool operator()(RandomGenerator& gen, F length) {
    return one_dis(gen) >= exp(-length * inv_scale);
  }

  template <typename RandomGenerator>
  _ F deposition_depth(RandomGenerator& gen) {
    auto r = one_dis(gen);
    return -log(r) * scale;
  }

  _ static F max_bias() { return static_cast<F>(0); }

 private:
  util::random::uniform_real_distribution<F> one_dis;
  F scale;
  F inv_scale;
};
}  // Barrel
}  // PET2D

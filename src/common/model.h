#pragma once

#include "util/random.h"
#include "util/cuda/compat.h"

/// Generic classes applicable to PET 2D & 3D and other tasks
namespace Common {

/// Model which always produces a decay

/// This model always produces a decay opposite to ScintillatorAccept.
template <typename FType> class AlwaysAccept {
 public:
  using F = FType;

  AlwaysAccept() {}
  template <class RNG> _ bool operator()(RNG&, F) { return true; }

  template <class RNG> _ F deposition_depth(RNG&) { return static_cast<F>(0); }
};

/// Model of scintilator acceptance

/// Represents model of scintilator where CDF of decay is given by:
/// \f[
///     F = 1-e^{-scale * length}
/// \f]
/// where for given length \f$l\f$ we get:
/// \f[
///     F(l) = \begin{cases}
///        0             & l = 0           \\&
///        1 - 1/sqrt(e) & l = 1/2 * scale \\&
///        1 - 1/e       & l = 1/scale     \\&
///        1 - 1/e^2     & l = 2/scale
///     \end{cases}
/// \f]
///
/// \sa AlwaysAccept

template <typename FType> class ScintillatorAccept {
 public:
  using F = FType;

  _ ScintillatorAccept(F scale)
      : one_dis(static_cast<F>(0), static_cast<F>(1)),
        scale(scale),
        inv_scale(static_cast<F>(1) / scale) {}

  template <class RNG> _ bool operator()(RNG& rng, F length) {
    return one_dis(rng) >= exp(-length * inv_scale);
  }

  template <class RNG> _ F deposition_depth(RNG& rng) {
    auto r = one_dis(rng);
    return -log(r) * scale;
  }

 private:
  util::random::uniform_real_distribution<F> one_dis;
  F scale;
  F inv_scale;
};
}  // Common

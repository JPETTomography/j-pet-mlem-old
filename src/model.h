#pragma once

#include "random.h"

template <typename F = double> class AlwaysAccept {
 public:
  AlwaysAccept() {}
  template <class RandomGenerator> bool operator()(RandomGenerator&, F) {
    return true;
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator&) {
    return 0.0;
  }
};

template <typename F = double> class ScintilatorAccept {
 public:
  ScintilatorAccept(F unit_prob)
      : one_dis_(0., 1.),
        unit_prob_(unit_prob),
        inv_unit_prob_(1.0 / unit_prob) {
  }

  template <class RandomGenerator>
  bool operator()(RandomGenerator& gen, F length) {
    return one_dis_(gen) >= exp(-length * unit_prob_);
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator& gen) {
    return -log(one_dis_(gen)) * inv_unit_prob_;
  }

 private:
  uniform_real_distribution<F> one_dis_;
  F unit_prob_;
  F inv_unit_prob_;
};

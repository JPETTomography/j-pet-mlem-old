#pragma once

#include "random.h"

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

template <typename FType = double> class ScintilatorAccept {
 public:
  typedef FType F;

  ScintilatorAccept(F unit_prob)
      : one_dis_(static_cast<F>(0), static_cast<F>(1)),
        unit_prob_(unit_prob),
        inv_unit_prob_(static_cast<F>(1) / unit_prob) {
  }

  template <class RandomGenerator>
  bool operator()(RandomGenerator& gen, F length) {
    return one_dis_(gen) >= exp(-length * unit_prob_);
  }

  template <typename RandomGenerator> F deposition_depth(RandomGenerator& gen) {
    return -log(one_dis_(gen)) * inv_unit_prob_;
  }

  static F max_bias() { return static_cast<F>(0); }

 private:
  uniform_real_distribution<F> one_dis_;
  F unit_prob_;
  F inv_unit_prob_;
};

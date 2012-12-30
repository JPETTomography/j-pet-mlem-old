#pragma once

#include "random.h"

template <typename F = double>
class always_accept {
public:
  always_accept() {}
  template <class RandomGenerator>
  bool operator () (RandomGenerator &, F) { return true; }

  template <typename RandomGenerator> 
  F deposition_depth(RandomGenerator &) {
    return 0.0;
  }
};

template <typename F = double>
class scintilator_accept {
public:
  scintilator_accept(F a_unit_prob): 
    one_dis(0., 1.),
    unit_prob(a_unit_prob),
    inv_unit_prob(1.0/unit_prob) {}

  template <class RandomGenerator>
  bool operator () (RandomGenerator &gen, F length) {
    return one_dis(gen) >= exp(-length * unit_prob);
  }

  template <typename RandomGenerator> 
  F deposition_depth(RandomGenerator &gen) {
    return -log(one_dis(gen))*inv_unit_prob;
  }

private:
  uniform_real_distribution<F> one_dis;
  F unit_prob;
  F inv_unit_prob;
};

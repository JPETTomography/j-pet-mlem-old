#pragma once

#include <random>

template <typename F = double>
class always_accept {
public:
  always_accept() {}
  template <class RandomGenerator>
  bool operator () (RandomGenerator &, F) { return true; }
};

template <typename F = double>
class scintilator_accept {
public:
  scintilator_accept(F a_unit_prob)
  : one_dis(0., 1.)
  , unit_prob(a_unit_prob)
  {}

  template <class RandomGenerator>
  bool operator () (RandomGenerator &gen, F length) {
    return one_dis(gen) >= exp(-length * unit_prob);
  }

private:
  std::uniform_real_distribution<F> one_dis;
  F unit_prob;
};

#pragma once

#include <random>

template <typename F = double>
class always_accept {
public:
  always_accept() {}
  bool operator () (F) { return true; }
};

template <class RandomGenerator, typename F = double>
class scintilator_accept {
public:
  scintilator_accept(RandomGenerator a_gen, F a_unit_prob)
  : gen(a_gen)
  , one_dis(0., 1.)
  , unit_prob(a_unit_prob)
  {}

  bool operator () (F length) {
    return one_dis(gen) >= exp(-length * unit_prob);
  }

private:
  RandomGenerator gen;
  std::uniform_real_distribution<F> one_dis;
  F unit_prob;
};

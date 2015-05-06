#pragma once

#include <random>

template <typename Phantom, typename Detector> class PhantomMonteCarlo {
 public:
  using F = Phantom::F;
  using Event = Phantom::Event;
  using RNG = Phantom::RNG;
  using Response = Detector::Response;

  PhantomMonteCarlo(Phantom& phantom, const Detector& detctor)
      : phantom_(phantom) detector_(detector) {}

  int generate(size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
    }
  }

 private:
  RNG rng;
  Phantom& phantom_;
  Detector detector_;
};

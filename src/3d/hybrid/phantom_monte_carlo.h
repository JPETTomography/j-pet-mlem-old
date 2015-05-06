#pragma once

#include <random>

template <typename Phantom, typename Detector> class PhantomMonteCarlo {
 public:
  using F = typename Phantom::F;
  using Event = typename Phantom::Event;
  using RNG = typename Phantom::RNG;
  using Response = typename Detector::Response;

  PhantomMonteCarlo(Phantom& phantom, const Detector& detector)
      : phantom_(phantom), detector_(detector) {}

  template<typename ModelType>
  int generate(RNG& rng, ModelType model, size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);

    }
    return 0;
  }

 private:
  Phantom& phantom_;
  Detector detector_;
};

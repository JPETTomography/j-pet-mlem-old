#pragma once

#include <ostream>
#include <iostream>
#include <random>
#include <functional>

namespace Common {

template <typename Phantom, typename Detector> class PhantomMonteCarlo {
 public:
  using F = typename Phantom::F;
  using Event = typename Phantom::Event;
  using RNG = typename Phantom::RNG;
  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  PhantomMonteCarlo(Phantom& phantom, const Detector& detector)
      : phantom_(phantom),
        detector_(detector),
        out_wo_error(std::cout),
        out_w_error(std::cout),
        out_exact_events(std::cout),
        out_full_response(std::cout) {}

  template <typename ModelType>
  int generate(RNG& rng, ModelType model, size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
      FullResponse full_response;

      if (detector_.exact_detect(rng, model, event, full_response) == 2) {
        out_full_response << full_response << std::endl;
        out_wo_error << detector_.noErrorResponse(full_response) << std::endl;
        out_w_error << detector_.errorResponse(rng, full_response) << std::endl;
        out_exact_events << event << std::endl;
      }
    }
    return 0;
  }

 private:
  Phantom& phantom_;
  Detector detector_;

 public:
  std::reference_wrapper<std::ostream> out_wo_error;
  std::reference_wrapper<std::ostream> out_w_error;
  std::reference_wrapper<std::ostream> out_exact_events;
  std::reference_wrapper<std::ostream> out_full_response;
};
}

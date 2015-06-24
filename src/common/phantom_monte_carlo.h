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

  static std::ofstream null_stream;

  PhantomMonteCarlo(Phantom& phantom, Detector& detector)
      : phantom_(phantom),
        detector_(detector),
        n_events_detected_(),
        out_wo_error(null_stream),
        out_w_error(null_stream),
        out_exact_events(null_stream),
        out_full_response(null_stream) {}

  template <typename ModelType>
  int generate(RNG& rng, ModelType model, size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
      FullResponse full_response;

      if (detector_.exact_detect(rng, model, event, full_response) == 2) {
        out_full_response << full_response << "\n";
        out_wo_error << detector_.response_wo_error(full_response) << "\n";
        out_w_error << detector_.response_w_error(rng, full_response) << "\n";
        out_exact_events << event << "\n";
        n_events_detected_++;
      }
    }
    return 0;
  }

  int n_events_detected() const { return n_events_detected_; }

 private:
  Phantom& phantom_;
  Detector& detector_;
  int n_events_detected_;

 public:
  std::reference_wrapper<std::ostream> out_wo_error;
  std::reference_wrapper<std::ostream> out_w_error;
  std::reference_wrapper<std::ostream> out_exact_events;
  std::reference_wrapper<std::ostream> out_full_response;
};

template <typename Phantom, typename Detector>
std::ofstream PhantomMonteCarlo<Phantom, Detector>::null_stream;
}

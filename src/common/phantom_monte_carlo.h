#pragma once

#include <random>
#include <functional>

namespace Common {

/// Executes Monte-Carlo for given \c Phantom and \c Detector classes

/// This is wrapper class that executes Monte-Carlo on abstract \c Phantom that
/// has to implement \c gen_event and \c Detector that has to implement \c
/// exact_detect method.
template <typename Phantom, typename Detector> class PhantomMonteCarlo {
 public:
  using F = typename Phantom::F;
  using Event = typename Phantom::Event;
  using RNG = typename Phantom::RNG;
  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  PhantomMonteCarlo(Phantom& phantom, Detector& detector)
      : phantom_(phantom), detector_(detector), n_events_detected_() {}

  template <class ModelClass, class EmitCallback, class DetectCallback>
  int operator()(RNG& rng,
                 ModelClass model,
                 size_t n_emisions,
                 EmitCallback emitted,
                 DetectCallback detected) {
    int n_events_detected = 0;
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
      emitted(event);
      FullResponse full_response;
      if (detector_.exact_detect(rng, model, event, full_response) == 2) {
        detected(event, full_response);
        n_events_detected++;
      }
    }
    n_events_detected_ += n_events_detected;
    return n_events_detected;
  }

  int n_events_detected() const { return n_events_detected_; }

 private:
  Phantom& phantom_;
  Detector& detector_;
  int n_events_detected_;
};
}

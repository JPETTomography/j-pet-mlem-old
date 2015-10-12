#pragma once

#include <random>
#include <functional>

namespace Common {

/// Executes Monte-Carlo for given \c Phantom and \c Detector classes
////
/// This is wrapper class that executes Monte-Carlo on abstract \c Phantom that
/// has to implement \c gen_event and \c Detector that has to implement \c
/// exact_detect method.
template <class PhantomType, class DetectorType, class ImageType>
class PhantomMonteCarlo {
 public:
  using Phantom = PhantomType;
  using Detector = DetectorType;
  using Image = ImageType;
  using F = typename Phantom::F;
  using Event = typename Phantom::Event;
  using RNG = typename Phantom::RNG;
  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  PhantomMonteCarlo(Phantom& phantom, Detector& detector)
      : phantom_(phantom),
        detector_(detector),
        n_events_emitted_(),
        n_events_detected_() {}

  template <class ModelClass,
            class EmitCallback,
            class DetectCallback,
            class Progress>
  int operator()(RNG& rng,
                 ModelClass model,
                 size_t n_emissions_or_detections,
                 EmitCallback emitted,
                 DetectCallback detected,
                 Progress progress,
                 bool only_detected = false) {
    size_t n_events_detected = 0;
    // (1) we run until we detected n_emissions_or_detections
    if (only_detected) {
      size_t n_events_emitted = 0;
      while (n_events_detected < n_emissions_or_detections) {
        progress(n_events_detected, false);
        auto event = phantom_.gen_event(rng);
        emitted(event);
        FullResponse full_response;
        if (detector_.exact_detect(rng, model, event, full_response) == 2) {
          detected(event, full_response);
          ++n_events_detected;
        }
        ++n_events_emitted;
        progress(n_events_detected, true);
      }
      n_events_emitted_ += n_events_emitted;
    }
    // (2) we run until we emit n_emissions_or_detections
    else {
      for (size_t n_events_emitted = 0;
           n_events_emitted < n_emissions_or_detections;
           ++n_events_emitted) {
        progress(n_events_emitted, false);
        auto event = phantom_.gen_event(rng);
        emitted(event);
        FullResponse full_response;
        if (detector_.exact_detect(rng, model, event, full_response) == 2) {
          detected(event, full_response);
          ++n_events_detected;
        }
        progress(n_events_emitted, true);
      }
      n_events_emitted_ += n_emissions_or_detections;
    }
    n_events_detected_ += n_events_detected;
    return n_events_detected;
  }

  int n_events_emitted() const { return n_events_emitted_; }
  int n_events_detected() const { return n_events_detected_; }

 private:
  Phantom& phantom_;
  Detector& detector_;
  size_t n_events_emitted_;
  size_t n_events_detected_;
};
}

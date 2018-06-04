#pragma once

#include <random>
#include <functional>

namespace Common {

/// Executes Monte-Carlo for given \c Phantom and \c Detector classes
////
/// This is wrapper class that executes Monte-Carlo on abstract \c Phantom that
/// has to implement \c gen_event and \c Detector that has to implement \c
/// exact_detect method.
template <class PhantomType, class DetectorType> class PhantomMonteCarlo {
 public:
  using Phantom = PhantomType;
  using Detector = DetectorType;
  using F = typename Phantom::F;
  using Event = typename Detector::Event;
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
#if _OPENMP
#pragma omp parallel shared(n_events_detected)
#endif
      {
#if _OPENMP
        RNG l_rng;
#pragma omp critical
        l_rng.seed(rng());
#else
#define l_rng rng
#endif
        while (n_events_detected < n_emissions_or_detections) {
          progress(n_events_detected, false);
          auto event = phantom_.gen_event(l_rng);
          emitted(event);
          FullResponse full_response;
          if (detector_.exact_detect(l_rng, model, event, full_response) == 2) {
            detected(event, full_response);
#if _OPENMP
#pragma omp atomic
#endif
            ++n_events_detected;
          }
          progress(n_events_detected, true);
        }
      }
      n_events_emitted_ += n_emissions_or_detections;
    }
    // (2) we run until we emit n_emissions_or_detections
    else {
#if _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
      for (
#if !_MSC_VER
          size_t n_events_emitted = 0;
          n_events_emitted < n_emissions_or_detections;
#else
          intptr_t n_events_emitted = 0;
          n_events_emitted < (intptr_t)n_emissions_or_detections;
#endif
          ++n_events_emitted) {
#if _OPENMP
        RNG l_rng;
#pragma omp critical
        l_rng.seed(rng());
#else
#define l_rng rng
#endif
        progress(n_events_emitted, false);
        auto event = phantom_.gen_event(l_rng);
        emitted(event);
        FullResponse full_response;
        if (detector_.exact_detect(l_rng, model, event, full_response) == 2) {
          detected(event, full_response);
#if _OPENMP
#pragma omp atomic
#endif
          ++n_events_detected;
        }
        progress(n_events_emitted, true);
      }
      n_events_emitted_ += n_emissions_or_detections;
    }
    n_events_detected_ += n_events_detected;
    return n_events_detected_;
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

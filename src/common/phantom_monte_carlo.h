#pragma once

#include <ostream>
#include <iostream>
#include <random>
#include <functional>

#include "util/null_stream.h"

namespace Common {

namespace {
static util::null_ostream null_ostream;
}

template <typename Phantom, typename Detector> class PhantomMonteCarlo {
 public:
  using F = typename Phantom::F;
  using Event = typename Phantom::Event;
  using RNG = typename Phantom::RNG;
  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  PhantomMonteCarlo(Phantom& phantom, Detector& detector)
      : phantom_(phantom),
        detector_(detector),
        n_events_detected_(),
        out_wo_error(null_ostream),
        out_w_error(null_ostream),
        out_exact_events(null_ostream),
        out_full_response(null_ostream) {}

  typename std::vector<FullResponse>::const_iterator begin() const {
    return responses_.begin();
  }

  typename std::vector<FullResponse>::const_iterator end() const {
    return responses_.end();
  }

  template <typename ModelType>
  int generate(RNG& rng, ModelType model, size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
      FullResponse full_response;

      if (detector_.exact_detect(rng, model, event, full_response) == 2) {
        responses_.push_back(full_response);
        out_exact_events << event << "\n";
        n_events_detected_++;
      }
    }
    return 0;
  }

  void write_out(RNG& rng) const {
    for (auto& full_response : responses_) {
      out_full_response << full_response << "\n";
      out_wo_error << detector_.response_wo_error(full_response) << "\n";
      out_w_error << detector_.response_w_error(rng, full_response) << "\n";
    }
  }

  int n_events_detected() const { return n_events_detected_; }

 private:
  Phantom& phantom_;
  Detector& detector_;
  int n_events_detected_;
  std::vector<FullResponse> responses_;

 public:
  std::reference_wrapper<std::ostream> out_wo_error;
  std::reference_wrapper<std::ostream> out_w_error;
  std::reference_wrapper<std::ostream> out_exact_events;
  std::reference_wrapper<std::ostream> out_full_response;
};
}

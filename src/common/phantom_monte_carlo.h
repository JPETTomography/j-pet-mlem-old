#pragma once
#include <ostream>

#include <random>

#include "3d/geometry/point.h"

namespace PET3D {

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
        no_error_stream(),
        error_stream(),
        exact_event_stream(),
        full_response_stream() {}

  void set_no_error_stream(std::ostream& os) { no_error_stream = &os; }
  void set_error_stream(std::ostream& os) { error_stream = &os; }
  void set_exact_event_stream(std::ostream& os) { exact_event_stream = &os; }
  void set_full_response_stream(std::ostream& os) {
    full_response_stream = &os;
  }

  template <typename ModelType>
  int generate(RNG& rng, ModelType model, size_t n_emisions) {
    for (size_t i = 0; i < n_emisions; ++i) {
      auto event = phantom_.gen_event(rng);
      FullResponse full_response;

      if (detector_.exact_detect(rng, model, event, full_response) == 2) {

        if (full_response_stream) {
          *full_response_stream << full_response << "\n";
        }

        if (no_error_stream) {
          Response response = detector_.noErrorResponse(full_response);
          *no_error_stream << response << "\n";
        }

        if (error_stream) {
          Response response = detector_.errorResponse(rng, full_response);
          *error_stream << response << "\n";
        }

        if (exact_event_stream)
          *exact_event_stream << event << "\n";
      }
    }
    return 0;
  }

 private:
  Phantom& phantom_;
  Detector detector_;
  std::ostream* no_error_stream;
  std::ostream* error_stream;
  std::ostream* exact_event_stream;
  std::ostream* full_response_stream;
};
}

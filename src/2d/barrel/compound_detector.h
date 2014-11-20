#pragma once

#include <vector>

#include "square_detector.h"
#include "circle_detector.h"
#include "util/svg_ostream.h"
#include "util/array.h"
#include "lor.h"
#ifndef MAX_DETECTORS
#define MAX_DETECTORS 256
#endif

#include <iostream>  // DEBUG

namespace PET2D {
namespace Barrel {

/// Detector made of several other detectors

/// No assumptions are made for how geometry of this detector looks like in
/// comparison to DetectorRing where are single detectors are placed on the
/// ring.
template <typename DetectorType = SquareDetector<double>,
          std::size_t MaxDetectors = MAX_DETECTORS,
          typename SType = int>
class CompoundDetector : public util::array<MaxDetectors, DetectorType> {
 public:
  using Detector = DetectorType;
  using S = SType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Pixel = PET2D::Pixel<S>;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;
  using Base = util::array<MaxDetectors, Detector>;
  using CircleDetector = Barrel::CircleDetector<F>;

  /// Tries to detect given event.

  /// \return number of coincidences (detector hits)
  /// \todo FIXME: Now returns first index to closest detector
  template <class RandomGenerator, class AcceptanceModel>
  short detect(RandomGenerator& gen,    ///< random number generator
               AcceptanceModel& model,  ///< acceptance model
               const Event& e,          ///< event to be detected
               LOR& lor,                ///<[out] lor of the event
               F& position              ///<[out] position of the event
               ) {
    (void)(gen);       // mark as used
    (void)(model);     // mark as used
    (void)(lor);       // mark as used
    (void)(position);  // mark as used
    int close_indices[MaxDetectors];
    int n_close_indices = 0;
    for (int i = 0; i < c_detectors.size(); ++i) {
      auto& circle = c_detectors[i];
      if (circle.intersects(e)) {
        close_indices[n_close_indices++] = i;
      }
    }
    std::sort(
        &close_indices[0], &close_indices[n_close_indices], [&](int a, int b) {
          return (e - c_detectors[close_indices[a]]).length2() <
                 (e - c_detectors[close_indices[b]]).length2();
        });
    return close_indices[0];
  }

  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          CompoundDetector& cd) {
    for (auto& detector : cd) {
      svg << detector;
    }

    return svg;
  }

  void push_back(const Detector& detector) {
    Base::push_back(detector);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  template <class... Args> void emplace_back(Args&&... args) {
    Base::emplace_back(std::forward<Args>(args)...);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  const CircleDetector& circumscribed(int i) { return c_detectors[i]; }

 private:
  util::array<MaxDetectors, CircleDetector> c_detectors;
};
}
}

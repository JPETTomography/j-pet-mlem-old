#pragma once

#include <vector>

#include "square_detector.h"
#include "util/svg_ostream.h"
#include "lor.h"

namespace PET2D {
namespace Barrel {

/// Detector made of several other detectors

/// No assumptions are made for how geometry of this detector looks like in
/// comparison to DetectorRing where are single detectors are placed on the
/// ring.
template <typename FType = double,
          typename SType = int,
          typename DetectorType = SquareDetector<FType>>
class CompoundDetector : std::vector<DetectorType> {
 public:
  typedef FType F;
  typedef SType S;
  typedef Barrel::LOR<S> LOR;
  typedef PET2D::Pixel<S> Pixel;
  typedef PET2D::Point<F> Point;
  typedef DetectorType Detector;
  typedef Barrel::Event<F> Event;

  /// Tries to detect given event.

  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  short detect(RandomGenerator& gen,    ///< random number generator
               AcceptanceModel& model,  ///< acceptance model
               const Event& e,          ///< event to be detected
               LOR& lor,                ///<[out] lor of the event
               F& position              ///<[out] position of the event
               ) {
    // FIXME: implement me!
    (void)(gen, model, e, lor, position);
    if (detector_centers.size() != this->size()) {
      calculate_detector_centers();
    }
    return 0;
  }

  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          CompoundDetector& cd) {
    for (auto& detector : cd) {
      svg << detector;
    }

    return svg;
  }

 private:
  std::vector<Point> detector_centers;

  void calculate_detector_centers() {
    detector_centers.resize(this->size());
    for (size_t i = 0; i < this->size(); ++i) {
      auto& center = detector_centers[i];
      center = (*this)[i]->center();
    }
  }
};
}
}

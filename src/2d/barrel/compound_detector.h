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
class CompoundDetector : public std::vector<DetectorType> {
  typedef FType F;
  typedef SType S;
  typedef Barrel::LOR<S> LOR;
  typedef PET2D::Pixel<S> Pixel;
  typedef PET2D::Point<F> Point;
  typedef DetectorType Detector;
  typedef Barrel::Event<F> Event;

  template <class RandomGenerator, class AcceptanceModel>
  short emit_event(RandomGenerator& gen,    ///< random number generator
                   AcceptanceModel& model,  ///< acceptance model
                   F rx,        ///< x coordinate of the emission point
                   F ry,        ///< y coordinate of the emission point
                   F angle,     ///< emission angle
                   LOR& lor,    ///<[out] lor of the event
                   F& position  ///<[out] position of the event
                   ) {
    // FIXME: implement me!
    (void)(gen, model, rx, ry, angle, lor, position);
    return 0;
  }

  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          CompoundDetector& cd) {
    for (auto& detector : cd) {
      svg << detector;
    }

    return svg;
  }
};
}
}

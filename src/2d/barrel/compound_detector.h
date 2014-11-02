#pragma once

#include <vector>

#include "square_detector.h"
#include "circle_detector.h"
#include "util/svg_ostream.h"
#include "lor.h"

namespace PET2D {
namespace Barrel {

/// Detector made of several other detectors

/// No assumptions are made for how geometry of this detector looks like in
/// comparison to DetectorRing where are single detectors are placed on the
/// ring.
template <typename DetectorType = SquareDetector<double>, typename SType = int>
class CompoundDetector : public std::vector<DetectorType> {
 public:
  typedef DetectorType Detector;
  typedef SType S;
  typedef typename Detector::F F;
  typedef Barrel::LOR<S> LOR;
  typedef PET2D::Pixel<S> Pixel;
  typedef PET2D::Point<F> Point;
  typedef Barrel::Event<F> Event;
  typedef std::vector<Detector> Base;
  typedef Barrel::CircleDetector<F> CircleDetector;

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
    return 0;
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
    c_detectors.push_back(this->back().circumscribe_center());
  }

  template <class... Args> void emplace_back(Args&&... args) {
    Base::emplace_back(std::forward<Args&&>(args)...);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  const CircleDetector& circumscribed(int i) { return c_detectors[i]; }

 private:
  std::vector<CircleDetector> c_detectors;
};
}
}

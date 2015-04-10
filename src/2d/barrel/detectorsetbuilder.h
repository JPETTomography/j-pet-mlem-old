#ifndef DETECTORSETBUILDER_H
#define DETECTORSETBUILDER_H

#include "2d/barrel/detector_set.h"

namespace PET2D {
namespace Barrel {

template <typename DetectorSetType> class DetectorSetBuilder {
 public:
  using F = typename DetectorSetType::F;
  using S = typename DetectorSetType::S;
  using Detector = typename DetectorSetType::Detector;
  using Vector = PET2D::Vector<F>;

  static DetectorSetType buildSingleRing(F radius,
                                         S n_detectors,
                                         F w_detector,
                                         F h_detector,
                                         F d_detector = 0) {

    if (n_detectors > static_cast<S>(DetectorSetType::MaxDetectors))
      throw("too many detectors");
    if (radius <= 0)
      throw("invalid radius");
    if (w_detector > 0 && h_detector == 0)
      h_detector = Detector::default_height_for_width(w_detector);
    // NOTE: detector may return 0 for default height, which means we need to
    // have height given explicitely
    if (w_detector <= 0 || h_detector <= 0)
      throw("invalid detector size");
    if (n_detectors % 4)
      throw("number of detectors must be multiple of 4");

    Detector detector_base(w_detector, h_detector, d_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base +=
        Vector(0, radius + (d_detector > 0 ? d_detector : h_detector) / 2);

    F outer_radius = radius + (d_detector > 0 ? d_detector : h_detector);
    if (d_detector == 0)
      outer_radius = detector_base.max_distance();

    DetectorSetType detector_set(radius, outer_radius);

    // produce detector ring rotating base detector n times
    for (auto n = 0; n < n_detectors; ++n) {
      auto detector = detector_base;
      detector.rotate(2 * F(M_PI) * n / n_detectors - F(M_PI_2));
      detector_set.push_back(detector);
    }

    return detector_set;
  }
};
}
}

#endif  // DETECTORSETBUILDER_H

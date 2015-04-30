#ifndef DETECTORSETBUILDER_H
#define DETECTORSETBUILDER_H

#include <functional>

#include "2d/barrel/detector_set.h"
#include "symmetry_descriptor.h"

namespace PET2D {
namespace Barrel {

template <typename DetectorSetType> class DetectorSetBuilder {
 public:
  using F = typename DetectorSetType::F;
  using S = typename DetectorSetType::S;
  using Detector = typename DetectorSetType::Detector;
  using Vector = PET2D::Vector<F>;

  static DetectorSetType buildSingleRing(F radius,
                                         int n_detectors,
                                         F w_detector,
                                         F h_detector,
                                         F d_detector = 0) {

    if (n_detectors > static_cast<S>(DetectorSetType::MaxDetectors))
      throw("buld single ring: too many detectors");
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

    auto symmetry_descriptor = new SymmetryDescriptor<S>(n_detectors, 8);
    for (S d = 0; d < n_detectors; ++d) {
      for (S s = 0; s < SymmetryDescriptor<S>::EIGHT; ++s) {
        // std::cerr << d << " " << s << " " << symmetric_detector(d, s) <<
        // "\n";
        symmetry_descriptor->set_symmetric_detector(
            d,
            s,
            symmetry_descriptor->ring_symmetric_detector(n_detectors, d, s));
      }
    }
    detector_set.symmetry_descriptor_ = symmetry_descriptor;

    return detector_set;
  }

  static DetectorSetType buildMultipleRings(
      const std::vector<F> radius,    ///< radiuses of ring
      const std::vector<F> rotation,  ///< rotation of each ring (0-1)
      std::vector<int> n_detectors,   ///< numbers of detectors on ring
      F w_detector,                   ///< width of single detector (along ring)
      F h_detector,                   ///< height/depth of single detector
                                      ///< (perpendicular to ring)
      F d_detector = 0                ///< diameter of circle single detector is
                                      ///< inscribed in
      ) {

    if (!radius.size())
      throw("must specify at least one radius");
    if (n_detectors.size() > radius.size())
      throw("number of numbers of detectors must be less or equal radiuses");

    bool symmetry_broken = false;
    for (int i = 0; i < radius.size(); ++i) {

      if (std::abs(rotation[i]) >= 1e-6 &&
          std::abs(rotation[i] - F(0.5)) > 1e-6) {
        std::cerr << "rotation = " << rotation[i] << "\n";
        symmetry_broken = true;
        break;
      }
    }

    if (symmetry_broken) {
      std::cerr << "Waring : symmetry is broken\n";
    }

    int total_n_detectors = 0;

    for (size_t i = 0; i < radius.size(); ++i) {
      if (!radius[i])
        break;
      if (!n_detectors[i])
        n_detectors[i] = n_detectors[i - 1];

      total_n_detectors += n_detectors[i];
    }

    auto symmetry_descriptor = new SymmetryDescriptor<S>(total_n_detectors, 8);

    S start_detector = 0;

    DetectorSetType detector_set = buildSingleRing(
        radius[0], n_detectors[0], w_detector, h_detector, d_detector);

    std::function<S(S, S)> symmetric_detector;

    if (std::abs(rotation[0]) < 1e-6) {
      symmetric_detector = [=](S d, S s) -> S {
        return symmetry_descriptor->ring_symmetric_detector(
            n_detectors[0], d, s);
      };
    } else if (std::abs(rotation[0] - 0.5) < 1e-6) {
      symmetric_detector = [=](S d, S s) -> S {
        return symmetry_descriptor->rotated_ring_symmetric_detector(
            n_detectors[0], d, s);
      };
    } else {
      symmetric_detector = [=](S, S) -> S { return 0; };
    }

    for (S d = 0; d < n_detectors[0]; ++d) {
      for (S s = 0; s < SymmetryDescriptor<S>::EIGHT; ++s) {
        // std::cerr << d << " " << s << " " << symmetric_detector(d, s) <<
        // "\n";
        symmetry_descriptor->set_symmetric_detector(
            d, s, symmetric_detector(d, s));
      }
    }
    start_detector += n_detectors[0];
    // Now create all following rings
    for (size_t i = 1; i < radius.size(); ++i) {
      if (!radius[i])
        break;
      if (!n_detectors[i])
        n_detectors[i] = n_detectors[i - 1];
      if (n_detectors[i] + detector_set.size() > DetectorSetType::MaxDetectors)
        throw("build multiple rings :too many detectors");

      DetectorSetLayout<Detector, DetectorSetType::MaxDetectors, S> ring =
          buildSingleRing(
              radius[i], n_detectors[i], w_detector, h_detector, d_detector);
      S detector_i = 0;

      if (std::abs(rotation[i]) < 1e-6) {
        symmetric_detector = [=](S d, S s) -> S {
          return symmetry_descriptor->ring_symmetric_detector(
              n_detectors[i], d, s);
        };
      } else if (std::abs(rotation[i] - 0.5) < 1e-6) {
        symmetric_detector = [=](S d, S s) -> S {
          return symmetry_descriptor->rotated_ring_symmetric_detector(
              n_detectors[i], d, s);
        };
      } else {
        symmetric_detector = [=](S, S) -> S { return 0; };
      }

      for (auto& detector : ring) {
        detector.rotate(2 * F(M_PI) * rotation[i] / ring.size());
        detector_set.push_back(detector);

        for (S s = 0; s < SymmetryDescriptor<S>::EIGHT; ++s) {
          //          std::cerr << start_detector << " " << detector_i << " " <<
          //          s << " "
          //                    << symmetric_detector(detector_i, s) << "\n";
          symmetry_descriptor->set_symmetric_detector(
              start_detector + detector_i,
              s,
              start_detector + symmetric_detector(detector_i, s));
        }
        detector_i++;
      }
      if (ring.outer_radius() > detector_set.outer_radius())
        detector_set.c_outer = ring.c_outer;
      start_detector += n_detectors[i];
    }

    detector_set.symmetry_descriptor_ = symmetry_descriptor;
    return detector_set;
  }
};
}
}

#endif  // DETECTORSETBUILDER_H

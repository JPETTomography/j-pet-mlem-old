#pragma once

#include <functional>
#include <limits>

#include "2d/barrel/generic_scanner.h"
#include "symmetry_descriptor.h"

namespace PET2D {
namespace Barrel {

/// Factory building scanners containing single or multiple rings of detector
template <class ScannerClass> class ScannerBuilder {
 public:
  using Scanner = ScannerClass;
  using F = typename Scanner::F;
  using S = typename Scanner::S;
  using Detector = typename Scanner::Detector;
  using Vector = PET2D::Vector<F>;

  static Scanner build_single_ring(F radius,
                                   int n_detectors,
                                   F w_detector,
                                   F h_detector,
                                   F d_detector = 0,
                                   F fov_radius = 0) {

    if (n_detectors > static_cast<S>(Scanner::MaxDetectors))
      throw("buld single ring: too many detectors");
    if (radius <= 0)
      throw("invalid radius");
    if (w_detector > 0 && h_detector == 0)
      h_detector = Detector::default_height_for_width(w_detector);
    // NOTE: detector may return 0 for default height, which means we need to
    // have height given explicitely
    if (w_detector <= 0 || h_detector <= 0)
      throw("invalid detector size");
    //    if (n_detectors % 4)
    //      throw("number of detectors must be multiple of 4");

    Detector detector_base(w_detector, h_detector, d_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base +=
        Vector(0, radius + (d_detector > 0 ? d_detector : h_detector) / 2);

    F outer_radius = radius + (d_detector > 0 ? d_detector : h_detector);
    if (d_detector == 0)
      outer_radius = detector_base.max_distance();

    Scanner detector_set(radius, outer_radius, fov_radius);

    // produce detector ring rotating base detector n times
    for (auto n = 0; n < n_detectors; ++n) {
      auto detector = detector_base;
      detector.rotate(2 * F(M_PI) * n / n_detectors - F(M_PI_2));
      detector_set.push_back(detector);
    }

    auto symmetry_descriptor = new SymmetryDescriptor<S>(n_detectors, 8);
    for (S d = 0; d < n_detectors; ++d) {
      for (S s = 0; s < SymmetryDescriptor<S>::EIGHT; ++s) {
        symmetry_descriptor->set_symmetric_detector(
            d,
            s,
            symmetry_descriptor->ring_symmetric_detector(n_detectors, d, s));
      }
    }
    detector_set.symmetry_descriptor_ = symmetry_descriptor;

    return detector_set;
  }

  static Scanner build_multiple_rings(
      const std::vector<F> radius,    ///< radiuses of ring
      const std::vector<F> rotation,  ///< rotation of each ring (0-1)
      std::vector<int> n_detectors,   ///< numbers of detectors on ring
      F w_detector,                   ///< width of single detector (along ring)
      F h_detector,                   ///< height/depth of single detector
                                      ///< (perpendicular to ring)
      F d_detector = 0,               ///< diameter of circle single detector is
                                      ///< inscribed in
      F fov_radius = 0                ///< field of view radius (0-automatic)
      ) {

    if (!radius.size())
      throw("must specify at least one radius");
    if (n_detectors.size() > radius.size())
      throw("number of numbers of detectors must be less or equal radiuses");

    for (unsigned int i = 0; i < radius.size(); ++i) {
      if (std::abs(rotation[i]) >= std::numeric_limits<F>::epsilon() &&
          std::abs(rotation[i] - F(0.5)) > std::numeric_limits<F>::epsilon()) {
        throw("broken symmetry");
      }
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

    Scanner detector_set = build_single_ring(radius[0],
                                             n_detectors[0],
                                             w_detector,
                                             h_detector,
                                             d_detector,
                                             fov_radius);

    std::function<S(S, S)> symmetric_detector;

    if (std::abs(rotation[0]) < std::numeric_limits<F>::epsilon()) {
      symmetric_detector = [=](S d, S s) -> S {
        return symmetry_descriptor->ring_symmetric_detector(
            n_detectors[0], d, s);
      };
    } else if (std::abs(rotation[0] - F(0.5)) <
               std::numeric_limits<F>::epsilon()) {
      symmetric_detector = [=](S d, S s) -> S {
        return symmetry_descriptor->rotated_ring_symmetric_detector(
            n_detectors[0], d, s);
      };
    } else {
      symmetric_detector = [=](S, S) -> S { return 0; };
    }

    for (S d = 0; d < n_detectors[0]; ++d) {
      for (S s = 0; s < SymmetryDescriptor<S>::EIGHT; ++s) {
#if DEBUG
        std::cerr << d << " " << s << " " << symmetric_detector(d, s) << "\n";
#endif
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
      if (n_detectors[i] + detector_set.size() > Scanner::MaxDetectors)
        throw("build multiple rings: too many detectors");

      DetectorSet<Detector, S, Scanner::MaxDetectors> ring = build_single_ring(
          radius[i], n_detectors[i], w_detector, h_detector, d_detector);
      S detector_i = 0;

      if (std::abs(rotation[i]) < std::numeric_limits<F>::epsilon()) {
        symmetric_detector = [=](S d, S s) -> S {
          return symmetry_descriptor->ring_symmetric_detector(
              n_detectors[i], d, s);
        };
      } else if (std::abs(rotation[i] - F(0.5)) <
                 std::numeric_limits<F>::epsilon()) {
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
#if DEBUG
          std::cerr << start_detector << " " << detector_i << " " << s << " "
                    << symmetric_detector(detector_i, s) << "\n";
#endif
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

}  // Barrel
}  // PET2D

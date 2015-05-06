#pragma once

#include "detector_set.h"
#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/geometry/circle.h"
#include "square_detector.h"
#include "circle_detector.h"
#include "lor.h"
#include "util/random.h"

namespace PET2D {
namespace Barrel {

/// Scanner made of one 2D ring of detectors

/// This is optimized GenericScanner using assumption \b all detectors lie on
/// ring, so some operations like possible secants can be done much quicker,
/// approximately 2x faster than using GenericScanner.
///
/// \image html detector_ring.pdf.png

template <typename DetectorType, std::size_t MaxDetectors, typename SType>
class RingScanner : public DetectorSet<DetectorType, MaxDetectors, SType> {
 public:
  using Base = DetectorSet<DetectorType, MaxDetectors, SType>;
  using S = SType;
  using Detector = DetectorType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Pixel = PET2D::Pixel<S>;
  using Circle = PET2D::Circle<F>;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;
  using Response = typename Base::Response;

#if DISABLE_CONSTRUCTOR
  Scanner(F radius,         ///< radius of ring
          S n_detectors,    ///< number of detectors on ring
          F w_detector,     ///< width of single detector (along ring)
          F h_detector,     ///< height/depth of single detector
                            ///< (perpendicular to ring)
          F d_detector = 0  ///< diameter of circle single detector is
                            ///< inscribed in
          )
      : Base(radius, n_detectors, w_detector, h_detector, d_detector) {
    if (n_detectors % 4)
      throw("number of detectors must be multiple of 4");
  }
#endif

  RingScanner(F radius = 1, F outer_radius = F(1.5))
      : Base(radius, outer_radius) {}

 private:
  template <class RandomGenerator, class AcceptanceModel>
  _ bool check_for_hits(RandomGenerator& gen,
                        AcceptanceModel& model,
                        S inner,
                        S outer,
                        Event e,
                        S& detector,
                        F& depth,
                        Point& p1,
                        Point& p2) const {

    const auto n_detectors = this->size();
    // tells in which direction we got shorter modulo distance
    S step = ((n_detectors + inner - outer) % n_detectors >
              (n_detectors + outer - inner) % n_detectors)
                 ? 1
                 : n_detectors - 1;
    S end = (outer + step) % n_detectors;
    for (auto i = inner; i != end; i = (i + step) % n_detectors) {
      auto intersections = (*this)[i].intersections(e);
      // check if we got 2 point intersection
      // then test the model against these points distance
      if (intersections.size() == 2) {
        auto deposition_depth = model.deposition_depth(gen);
#if DEBUG
        std::cerr << "dep " << deposition_depth << " "
                  << (intersections[1] - intersections[0]).length()
                  << std::endl;
#endif
        if (deposition_depth < (intersections[1] - intersections[0]).length()) {
          detector = i;
          depth = deposition_depth;
          p1 = intersections[0];
          p2 = intersections[1];
          return true;
        }
      }
    }
    return false;
  }

 public:
  /// Tries to detect given event.

  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 Response& response       ///< scanner response (LOR+length)
                 ) const {

    const auto n_detectors = this->size();
#if !_MSC_VER
    const auto& c_inner = this->c_inner;
    const auto& c_outer = this->c_outer;
#endif
    auto inner_secant = c_inner.secant(e);
    auto outer_secant = c_outer.secant(e);

    if (inner_secant.size() != 2 || outer_secant.size() != 2)
      return 0;

    auto i_inner = c_inner.section(c_inner.angle(inner_secant[0]), n_detectors);
    auto i_outer = c_outer.section(c_inner.angle(outer_secant[0]), n_detectors);
    S detector1;
    F depth1;

    Point d1_p1, d1_p2;
    if (!check_for_hits(
            gen, model, i_inner, i_outer, e, detector1, depth1, d1_p1, d1_p2))
      return 0;

    i_inner = c_inner.section(c_inner.angle(inner_secant[1]), n_detectors);
    i_outer = c_outer.section(c_inner.angle(outer_secant[1]), n_detectors);
    S detector2;
    F depth2;
    Point d2_p1, d2_p2;
    if (!check_for_hits(
            gen, model, i_inner, i_outer, e, detector2, depth2, d2_p1, d2_p2))
      return 0;

    response.lor = LOR(detector1, detector2);

#if !__CUDACC__
    // FIXME: consider removing this check
    if (response.lor.first == response.lor.second)
      throw("invalid LOR");
#endif

    Point origin(e.x, e.y);
    F length1 = origin.nearest_distance(d1_p1, d1_p2) + depth1;
    F length2 = origin.nearest_distance(d2_p1, d2_p2) + depth2;

    if (detector1 > detector2) {
      response.dl = length1 - length2;
    } else {
      response.dl = length2 - length1;
    }

    return 2;
  }
};
}  // Barrel
}  // PET2D

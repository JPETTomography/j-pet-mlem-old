#pragma once

#include "square_detector.h"
#include "circle_detector.h"
#include "util/array.h"
#include "util/sort.h"
#include "lor.h"
#include "detector_set.h"

#if !__CUDACC__
#include "util/svg_ostream.h"
#include <vector>  // multi-ring detector construction
#endif

/// Two-dimensional PET
namespace PET2D {
/// Two-dimensional PET barrel
namespace Barrel {

/// Scanner made of several detectors

/// Represents scanner (PET) made of several detectors (scintillators)
/// using any geometry, particuallary it may be one or several rings of
/// detectors, or even detectors organized into other shapes.
///
/// No assumptions are made for how geometry of this scanner looks like in
/// comparison to ScannerRing where are single detectors are placed on the
/// ring.
///
/// \image html config_4x48.pdf.png
///
/// This class is a template that accepts custom DetectorType which can be any
/// shape of:
/// - SquareDetector
/// - CircleDetector
/// - TriangleDetector
/// - PolygonalDetector

template <typename DetectorType, std::size_t MaxDetectors, typename SType>
class GenericScanner : public DetectorSet<DetectorType, MaxDetectors, SType> {
 public:
  using Detector = DetectorType;
  using S = SType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Circle = PET2D::Circle<F>;
  using Event = Barrel::Event<F>;
  using Base = DetectorSet<DetectorType, MaxDetectors, SType>;
  using CircleDetector = Barrel::CircleDetector<F>;
  using Indices = util::array<MaxDetectors, S>;
  using Response = typename Base::Response;

  /// Makes an empty detector set.
  GenericScanner(F radius = 1, F outer_radius = F(1.5))
      : Base(radius, outer_radius) {}

  enum class TestCase {
    TEST_8_SQUARE_DETECTORS,
  };

  /// Makes new detector using hardcoded test case
  GenericScanner(TestCase test_case,  ///< test case
                 F radius,            ///< radius of ring
                 F w_detector,        ///< width of single detector (along ring)
                 F h_detector,        ///< height/depth of single detector
                                      ///< (perpendicular to ring)
                 F d_detector = 0     ///< diameter of circle single detector is
                                      ///< inscribed in
                 )
      : Base(test_case, radius, w_detector, h_detector, d_detector) {}

  /// \brief Tries to detect given event.
  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 Response& response       ///< scanner response (LOR+length)
                 ) const {
    Indices left, right;
    close_indices(e, left, right);
    S detector1, detector2;
    F depth1, depth2;
    Point d1_p1, d1_p2, d2_p1, d2_p2;
    if (!check_for_hits(gen, model, left, e, detector1, depth1, d1_p1, d1_p2) ||
        !check_for_hits(gen, model, right, e, detector2, depth2, d2_p1, d2_p2))
      return 0;

    response.lor = LOR(detector1, detector2);

    Point origin(e.center.x, e.center.y);
    F length1 = origin.nearest_distance(d1_p1, d1_p2) + depth1;
    F length2 = origin.nearest_distance(d2_p1, d2_p2) + depth2;

    if (detector1 > detector2) {
      response.dl = length1 - length2;
    } else {
      response.dl = length2 - length1;
    }

    return 2;
  }

  /// Produce indices of detectors close to given event
  _ void close_indices(const Event& e,     ///< event to be detected
                       Indices& negative,  ///<[out] indices on one side
                       Indices& positive   ///<[out] indices other side
                       ) const {
    F distances[MaxDetectors];
    auto pe = e.perpendicular();
    // select only these crossing circle circumscribed on detector
    for (int i = 0; i < static_cast<int>(Base::c_detectors.size()); ++i) {
      auto& circle = Base::c_detectors[i];
      if (circle.intersects(e)) {
        auto distance = pe(circle.center);
        distances[i] = distance;
        if (distance < 0) {
          negative.emplace_back(i);
        } else {
          positive.emplace_back(i);
        }
      }
    }
    // sort them so the closest go first
    // (1) negative distances (one side)
    util::heap_sort(negative.begin(),
                    negative.end(),
                    [&](S a, S b) { return distances[a] > distances[b]; });
    // (2) positive distances (other side)
    util::heap_sort(positive.begin(),
                    positive.end(),
                    [&](S a, S b) { return distances[a] < distances[b]; });
  }

  _ bool did_intersect(Event e, S detector, Point& p1, Point& p2) const {

    auto intersections = (*this)[detector].intersections(e);
    // check if we got 2 point intersection
    // then test the model against these points distance
    if (intersections.size() == 2) {

      p1 = intersections[0];
      p2 = intersections[1];
      return true;
    }

    return false;
  }

 private:
  template <class RandomGenerator, class AcceptanceModel>
  _ bool did_deposit(RandomGenerator& gen,
                     AcceptanceModel& model,
                     const Point& p1,
                     const Point& p2,
                     F& depth) const {
    depth = model.deposition_depth(gen);
#if DEBUG
    std::cerr << "dep " << depth << " " << (p2 - p1).length() << std::endl;
#endif
    if (depth < (p2 - p1).length()) {
      return true;
    }
    return false;
  }

  template <class RandomGenerator, class AcceptanceModel>
  _ bool check_for_hits(RandomGenerator& gen,
                        AcceptanceModel& model,
                        const Indices& indices,
                        Event e,
                        S& detector,
                        F& depth,
                        Point& p1,
                        Point& p2) const {

    for (auto i : indices) {
      if (did_intersect(e, i, p1, p2)) {
        if (did_deposit(gen, model, p1, p2, depth)) {
          detector = i;
          return true;
        }
      }
    }
    return false;
  }
};
}
}

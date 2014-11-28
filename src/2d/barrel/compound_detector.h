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
  using Indices = util::array<MaxDetectors, S>;
  using SideIndices = std::pair<Indices, Indices>;

  /// Tries to detect given event.

  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 LOR& lor,                ///<[out] lor of the event
                 F& position              ///<[out] position of the event
                 ) {
    auto indices = close_indices(e);
    S detector1, detector2;
    F depth1, depth2;
    Point d1_p1, d1_p2, d2_p1, d2_p2;
    if (!check_for_hits(
            gen, model, indices.first, e, detector1, depth1, d1_p1, d1_p2) ||
        !check_for_hits(
            gen, model, indices.second, e, detector2, depth2, d2_p1, d2_p2))
      return 0;

    lor = LOR(detector1, detector2);

    Point origin(e.x, e.y);
    F length1 = origin.nearest_distance(d1_p1, d1_p2) + depth1;
    F length2 = origin.nearest_distance(d2_p1, d2_p2) + depth2;

    if (detector1 > detector2) {
      position = length1 - length2;
    } else {
      position = length2 - length1;
    }

    return 2;
  }

  /// \return indices of detectors close to given event
  SideIndices close_indices(const Event& e  ///< event to be detected
                            ) {
    SideIndices indices;
    S distances[MaxDetectors];
    auto pe = e.perpendicular();
    // select only these crossing circle circumscribed on detector
    for (int i = 0; i < c_detectors.size(); ++i) {
      auto& circle = c_detectors[i];
      if (circle.intersects(e)) {
        auto distance = pe(circle);
        distances[i] = distance;
        if (distance < 0)
          indices.first.emplace_back(i);
        else
          indices.second.emplace_back(i);
      }
    }
    // sort them so the closest go first
    std::sort(indices.first.begin(), indices.first.end(), [&](S a, S b) {
      return distances[a] > distances[b];
    });
    std::sort(indices.second.begin(), indices.second.end(), [&](S a, S b) {
      return distances[a] < distances[b];
    });
    return indices;
  }

 private:
  template <class RandomGenerator, class AcceptanceModel>
  _ bool check_for_hits(RandomGenerator& gen,
                        AcceptanceModel& model,
                        const Indices& indices,
                        Event e,
                        S& detector,
                        F& depth,
                        Point& p1,
                        Point& p2) {

    for (auto i : indices) {
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

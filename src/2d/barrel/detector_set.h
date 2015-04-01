#pragma once

#include "square_detector.h"
#include "circle_detector.h"
#include "util/array.h"
#include "util/sort.h"
#include "lor.h"
#if !__CUDACC__
#include "util/svg_ostream.h"
#include <vector>  // multi-ring detector construction
#endif
#ifndef MAX_DETECTORS
#define MAX_DETECTORS 256
#endif

/// Two-dimensional PET
namespace PET2D {
/// Two-dimensional PET barrel
namespace Barrel {

/// Detector made of several other detectors

/// Represents detector compound made of several detectors (scintillators)
/// using any geometry, particuallary it may be one or several rings of
/// detectors, or even detectors organized into other shapes.
///
/// No assumptions are made for how geometry of this detector looks like in
/// comparison to DetectorRing where are single detectors are placed on the
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

template <typename DetectorType = SquareDetector<double>,
          std::size_t MaxDetectors = MAX_DETECTORS,
          typename SType = int>
class DetectorSet : public util::array<MaxDetectors, DetectorType> {
 public:
  using Detector = DetectorType;
  using S = SType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Point = PET2D::Point<F>;
  using Vector = PET2D::Vector<F>;
  using Circle = PET2D::Circle<F>;
  using Event = Barrel::Event<F>;
  using Base = util::array<MaxDetectors, Detector>;
  using CircleDetector = Barrel::CircleDetector<F>;
  using Indices = util::array<MaxDetectors, S>;

  /// Makes an empty detector set.
  DetectorSet(F radius = 1, F h_detector = 1)
      : Base(),
        fov_radius(radius / M_SQRT2),
        c_inner(radius),
        c_outer(radius + h_detector) {}

  /// Makes new detector set with detectors placed on the ring of given radius.
  DetectorSet(F radius,         ///< radius of ring
              S n_detectors,    ///< number of detectors on ring
              F w_detector,     ///< width of single detector (along ring)
              F h_detector,     ///< height/depth of single detector
                                ///< (perpendicular to ring)
              F d_detector = 0  ///< diameter of circle single detector is
                                ///< inscribed in
              )
      : Base(),
        fov_radius(radius / M_SQRT2),
        c_inner(radius),
        c_outer(radius + (d_detector > 0 ? d_detector : h_detector)) {
    if (n_detectors > static_cast<S>(MaxDetectors))
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

    // fix up outer circle
    if (d_detector == 0) {
      c_outer = Circle(detector_base.max_distance());
    }

    // produce detector ring rotating base detector n times
    for (auto n = 0; n < n_detectors; ++n) {
      auto detector = detector_base;
      detector.rotate(2 * F(M_PI) * n / n_detectors - F(M_PI_2));
      this->push_back(detector);
    }
  }

#if !__CUDACC__
  /// Makes new detector set with several rings.
  DetectorSet(const std::vector<F> radius,    ///< radiuses of ring
              const std::vector<F> rotation,  ///< rotation of each ring (0-1)
              std::vector<S> n_detectors,     ///< numbers of detectors on ring
              F w_detector,     ///< width of single detector (along ring)
              F h_detector,     ///< height/depth of single detector
                                ///< (perpendicular to ring)
              F d_detector = 0  ///< diameter of circle single detector is
                                ///< inscribed in
              )
      : DetectorSet(radius[0],
                    n_detectors[0],
                    w_detector,
                    h_detector,
                    d_detector) {
    if (!radius.size())
      throw("must specify at least one radius");
    if (n_detectors.size() > radius.size())
      throw("number of numbers of detectors must be less or equal radiuses");

    // Now create all following rings
    for (size_t i = 1; i < radius.size(); ++i) {
      if (!radius[i])
        break;
      if (!n_detectors[i])
        n_detectors[i] = n_detectors[i - 1];
      if (n_detectors[i] + this->size() > MaxDetectors)
        throw("too many detectors");
      DetectorSet ring(radius[i], n_detectors[i], w_detector, d_detector);
      for (auto& detector : ring) {
        detector.rotate(2 * F(M_PI) * rotation[i] / ring.size());
        this->push_back(detector);
      }
      if (ring.outer_radius() > outer_radius())
        c_outer = ring.c_outer;
    }
  }
#endif

  enum class TestCase {
    TEST_8_SQUARE_DETECTORS,
  };

  /// Makes new detector using hardcoded test case
  DetectorSet(TestCase test_case,  ///< test case
              F radius,            ///< radius of ring
              F w_detector,        ///< width of single detector (along ring)
              F h_detector,        ///< height/depth of single detector
                                   ///< (perpendicular to ring)
              F d_detector = 0     ///< diameter of circle single detector is
                                   ///< inscribed in
              )
      : Base(),
        fov_radius(radius / M_SQRT2),
        c_inner(radius),
        c_outer(radius + (d_detector > 0 ? d_detector : h_detector)) {

    Detector detector_base(w_detector, h_detector, d_detector);

    switch (test_case) {
      case TestCase::TEST_8_SQUARE_DETECTORS:
        this->push_back(detector_base + Point(-radius, -radius));
        this->push_back(detector_base + Point(-radius, 0));
        this->push_back(detector_base + Point(-radius, radius));
        this->push_back(detector_base + Point(radius, -radius));
        this->push_back(detector_base + Point(radius, 0));
        this->push_back(detector_base + Point(radius, radius));
        this->push_back(detector_base + Point(0, -radius));
        this->push_back(detector_base + Point(0, radius));
        break;
      default:
        throw("unknown test case");
    }
  }

  F radius() const { return c_inner.radius; }
  F outer_radius() const { return c_outer.radius; }
  F max_dl(F max_bias_size) const { return 2 * c_outer.radius + max_bias_size; }

  /// \brief Tries to detect given event.
  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 LOR& lor,                ///<[out] lor of the event
                 F& position              ///<[out] position of the event
                 ) const {
    Indices left, right;
    close_indices(e, left, right);
    S detector1, detector2;
    F depth1, depth2;
    Point d1_p1, d1_p2, d2_p1, d2_p2;
    if (!check_for_hits(gen, model, left, e, detector1, depth1, d1_p1, d1_p2) ||
        !check_for_hits(gen, model, right, e, detector2, depth2, d2_p1, d2_p2))
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

  /// Quantizes position across lor
  _ static S quantize_tof_position(F position,    ///< position across lor
                                   F step_size,   ///< step size
                                   S n_positions  ///< number of positions
                                   ) {
    // number of positions if always even, lower half are negative positions
    // where 0 means position closests to detector with higher index
    // maximum means position closests to detector with lower index
    if (position < 0)
      return n_positions / 2 - 1 -
             static_cast<S>(compat::floor(-position / step_size));
    else
      return static_cast<S>(compat::floor(position / step_size)) +
             n_positions / 2;
  }

  /// Returns number of position steps (indexes)
  S n_tof_positions(F step_size,     ///< step size
                    F max_bias_size  ///< possible bias (fuzz) maximum size
                    ) const {
    // since position needs to be symmetric against (0,0) number must be even
    return (static_cast<S>(ceil(2 * max_dl(max_bias_size) / step_size)) + 1) /
           2 * 2;
  }

  /// Produce indices of detectors close to given event
  _ void close_indices(const Event& e,     ///< event to be detected
                       Indices& negative,  ///<[out] indices on one side
                       Indices& positive   ///<[out] indices other side
                       ) const {
    F distances[MaxDetectors];
    auto pe = e.perpendicular();
    // select only these crossing circle circumscribed on detector
    for (int i = 0; i < static_cast<int>(c_detectors.size()); ++i) {
      auto& circle = c_detectors[i];
      if (circle.intersects(e)) {
        auto distance = pe(circle);
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

 private:
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

 public:
#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          DetectorSet& cd) {
    svg << cd.c_outer;
    svg << cd.c_inner;

    svg << "<g id=\"photomultipiers\">" << std::endl;
    for (auto& detector : cd.c_detectors) {
      svg << detector;
    }
    svg << "</g>" << std::endl;

    svg << "<g id=\"scintillators\">" << std::endl;
    for (auto& detector : cd) {
      svg << detector;
    }
    svg << "</g>" << std::endl;

    return svg;
  }
#endif

  void push_back(const Detector& detector) {
    Base::push_back(detector);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  void push_back(Detector&& detector) {
    Base::push_back(std::forward<Detector>(detector));
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  template <class... Args> void emplace_back(Args&&... args) {
    Base::emplace_back(std::forward<Args>(args)...);
    c_detectors.push_back(this->back().circumscribe_circle());
  }

  const CircleDetector& circumscribed(int i) const { return c_detectors[i]; }

  const F fov_radius;

 protected:
  util::array<MaxDetectors, CircleDetector> c_detectors;
  Circle c_inner;
  Circle c_outer;
};
}
}

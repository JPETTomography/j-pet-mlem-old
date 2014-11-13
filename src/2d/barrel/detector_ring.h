#pragma once

#include "util/array.h"
#include "2d/geometry/point.h"
#include "2d/geometry/pixel.h"
#include "2d/geometry/circle.h"
#include "square_detector.h"
#include "circle_detector.h"
#include "lor.h"
#include "util/random.h"
#if !__CUDACC__
#include "util/svg_ostream.h"
#endif

namespace PET2D {
namespace Barrel {

/// Detector made of 2D ring of single detectors

/// This is optimized CompoundDetector using assumption all detectors lie on
/// ring, so some operations like possible secants can be done much quicker.
template <typename DetectorType = SquareDetector<double>,
          typename SType = int,
          std::size_t MaxDetectors = 512>
class DetectorRing : public util::array<MaxDetectors, DetectorType> {
 public:
  using S = SType;
  using Detector = DetectorType;
  using F = typename Detector::F;
  using LOR = Barrel::LOR<S>;
  using Pixel = PET2D::Pixel<S>;
  using Circle = PET2D::Circle<F>;
  using Point = PET2D::Point<F>;
  using Event = Barrel::Event<F>;

  DetectorRing(S n_detectors,      ///< number of detectors on ring
               F radius,           ///< radius of ring
               F w_detector,       ///< width of single detector (along ring)
               F h_detector,       ///< height/depth of single detector
                                   ///< (perpendicular to ring)
               F d_detector = F()  ///< diameter of circle single detector is
                                   ///< inscribed in
               )
      : c_inner(radius),
        c_outer(radius + (d_detector > F() ? d_detector : h_detector)),
        n_detectors(n_detectors),
        n_lors(n_detectors * (n_detectors + 1) / 2),
        radius_diff(c_outer.radius - c_inner.radius) {
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

    fov_radius_ = radius / M_SQRT2;

    Detector detector_base(w_detector, h_detector, d_detector);
    auto r_detector = d_detector / 2;
    CircleDetector<F> circle_detector_base(r_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base +=
        Point(0, radius + (d_detector > F() ? d_detector : h_detector) / 2);
    circle_detector_base += Point(0, radius + r_detector);

    // fix up outer circle
    if (d_detector == F()) {
      c_outer = Circle(detector_base.max_distance());
    }

    // produce detector ring rotating base detector n times
    for (auto n = 0; n < n_detectors; ++n) {
      auto detector = detector_base;
      detector.rotate(2 * F(M_PI) * n / n_detectors - F(M_PI_2));
      this->push_back(detector);
      if (d_detector > F()) {
        auto circle_detector = circle_detector_base;
        circle_detector.rotate(2 * F(M_PI) * n / n_detectors - F(M_PI_2));
        c_detectors.push_back(circle_detector);
      }
    }
  }

  F radius() const { return c_inner.radius; }
  F outer_radius() const { return c_outer.radius; }
  S lors() const { return n_lors; }
  S detectors() const { return n_detectors; }
  F fov_radius() const { return fov_radius_; }

  F max_dl(F max_bias_size) const {
    return 2.0 * c_outer.radius + max_bias_size;
  }

  /// Quantizes position across lor
  _ S quantize_position(F position,    ///< position across lor
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
  S n_positions(F step_size,     ///< step size
                F max_bias_size  ///< possible bias (fuzz) maximum size
                ) {
    // since position needs to be symmetric against (0,0) number must be even
    return (static_cast<S>(ceil(2.0 * max_dl(max_bias_size) / step_size)) + 1) /
           2 * 2;
  }

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
                        Point& p2) {

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

  template <class RandomGenerator, class AcceptanceModel>
  _ bool check_for_hits(RandomGenerator& gen,
                        AcceptanceModel& model,
                        S inner,
                        S outer,
                        Event e,
                        S& detector,
                        F& depth) {
    Point p1, p2;
    return check_for_hits(gen, model, inner, outer, e, detector, depth, p1, p2);
  }

 public:
  /// Tries to detect given event.

  /// \return number of coincidences (detector hits)
  template <class RandomGenerator, class AcceptanceModel>
  _ short detect(RandomGenerator& gen,    ///< random number generator
                 AcceptanceModel& model,  ///< acceptance model
                 const Event& e,          ///< event to be detected
                 LOR& lor,                ///<[out] lor of the event
                 F& position              ///<[out] position of the event
                 ) {

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

    lor = LOR(detector1, detector2);

#if !__CUDACC__
    // FIXME: consider removing this check
    if (lor.first == lor.second)
      throw("invalid LOR");
#endif

    Point origin(e.x, e.y);
    F length1 = origin.nearest_distance(d1_p1, d1_p2) + depth1;
    F length2 = origin.nearest_distance(d2_p1, d2_p2) + depth2;

    if (detector1 > detector2) {
      position = length1 - length2;
    } else {
      position = length2 - length1;
    }

#ifdef GPU_TOF_TEST
    printf("%d %d %f %f %f %f %f %f %f %f %f\n",
           detector1,
           detector2,
           d1_p1.x,
           d1_p1.y,
           d1_p2.x,
           d1_p2.y,
           d2_p1.x,
           d2_p1.y,
           d2_p2.x,
           d2_p2.y,
           angle);
    printf("Length1: %f Length2: $%f\n",
           origin.nearest_distance(d1_p1, d1_p2),
           origin.nearest_distance(d2_p1, d2_p2));

#endif
#if DEBUG
    std::cerr << "position " << position << std::endl;
#endif
    return 2;
  }

#if !__CUDACC__
  friend util::svg_ostream<F>& operator<<(util::svg_ostream<F>& svg,
                                          DetectorRing& dr) {
    svg << dr.c_outer;
    svg << dr.c_inner;

    svg << "<g id=\"photomultipiers\">" << std::endl;
    for (auto& detector : dr.c_detectors) {
      svg << detector;
    }
    svg << "</g>" << std::endl;

    svg << "<g id=\"scintillators\">" << std::endl;
    for (auto& detector : dr) {
      svg << detector;
    }
    svg << "</g>" << std::endl;

    return svg;
  }
#endif

 private:
  Circle c_inner;
  Circle c_outer;

 public:
  const S n_detectors;
  const S n_lors;

 private:
  util::array<MaxDetectors, CircleDetector<F>> c_detectors;
  F fov_radius_;
  F radius_diff;
};
}  // Barrel
}  // PET2D

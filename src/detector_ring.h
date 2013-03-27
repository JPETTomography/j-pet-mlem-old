#pragma once

#include <map>

#include "random.h"
#include "detector.h"
#include "circle.h"
#include "svg_ostream.h"
#include "point.h"
#include "pixel.h"
#include "lor.h"

/// Provides model for 2D ring of detectors
template <typename FType = double, typename SType = int>
class DetectorRing : public std::vector<Detector<FType>> {
 public:
  typedef FType F;
  typedef SType S;
  typedef ::LOR<S> LOR;
  typedef ::Pixel<S> Pixel;
  typedef ::Circle<F> Circle;
  typedef ::Point<F> Point;
  typedef ::Detector<F> Detector;
  typedef ::Event<F> Event;

  /// @param n_detectors number of detectors on ring
  /// @param radius      radius of ring
  /// @param w_detector  width of single detector (along ring)
  /// @param h_detector  height/depth of single detector
  ///                    (perpendicular to ring)
  DetectorRing(S n_detectors, F radius, F w_detector, F h_detector)
      : c_inner_(radius),
        c_outer_(radius + h_detector),
        n_detectors_(n_detectors),
        n_lors_(n_detectors * (n_detectors + 1) / 2),
        radius_diff_(c_outer_.radius() - c_inner_.radius()) {
    if (radius <= 0.)
      throw("invalid radius");
    if (w_detector <= 0. || h_detector <= 0.)
      throw("invalid detector size");
    if (n_detectors_ % 4)
      throw("number of detectors must be multiple of 4");

    fov_radius_ = radius / M_SQRT2;

    Detector detector_base(h_detector, w_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base += Point(radius + h_detector / 2, 0.);

    // produce detector ring rotating base detector n times
    for (auto n = 0; n < n_detectors_; ++n) {
      this->push_back(detector_base.rotated(2. * M_PI * n / n_detectors_));
    }
  }

  F radius() const { return c_inner_.radius(); }
  S lors() const { return n_lors_; }
  S detectors() const { return n_detectors_; }
  F fov_radius() const { return fov_radius_; }

  Pixel pixel(F x, F y, F pixel_size) {
    F rx = x + fov_radius();
    F ry = y + fov_radius();
    return Pixel(static_cast<S>(floor(rx / pixel_size)),
                 static_cast<S>(floor(ry / pixel_size)));
  }

  F max_dl(F max_bias_size) const {
    return 2.0 * c_outer_.radius() + max_bias_size;
  }

  /// Quantizes position with given:
  /// @param step_size      step size
  /// @param max_bias_size  possible bias (fuzz) maximum size
  S quantize_position(F position, F step_size, F max_bias_size) {
    // FIXME: rounding?
    return static_cast<S>(
        floor((position + max_dl(max_bias_size)) / step_size));
  }

  /// Returns number of position steps (indexes) for:
  /// @param step_size      step size
  /// @param max_bias_size  possible bias (fuzz) maximum size
  S n_positions(F step_size, F max_bias_size) {
    return static_cast<S>(ceil(2.0 * max_dl(max_bias_size) / step_size));
  }

  template <class RandomGenerator, class AcceptanceModel>
  bool check_for_hits(RandomGenerator& gen,
                      AcceptanceModel& model,
                      S inner,
                      S outer,
                      Event e,
                      S& detector,
                      F& depth,
                      Point& p1,
                      Point& p2) {

    // tells in which direction we got shorter modulo distance
    S step = ((n_detectors_ + inner - outer) % n_detectors_ >
              (n_detectors_ + outer - inner) % n_detectors_)
             ? 1
             : n_detectors_ - 1;
    S end = (outer + step) % n_detectors_;
    for (auto i = inner; i != end; i = (i + step) % n_detectors_) {
      auto points = (*this)[i].intersections(e);
      // check if we got 2 point intersection
      // then test the model against these points distance
      if (points.size() == 2) {
        auto deposition_depth = model.deposition_depth(gen);
#if DEBUG
        std::cerr << "dep " << deposition_depth << " "
                  << (points[1] - points[0]).length() << std::endl;
#endif
        if (deposition_depth < (points[1] - points[0]).length()) {
          detector = i;
          depth = deposition_depth;
          p1 = points[0];
          p2 = points[1];
          return true;
        }
      }
    }
    return false;
  }

  template <class RandomGenerator, class AcceptanceModel>
  bool check_for_hits(RandomGenerator& gen,
                      AcceptanceModel& model,
                      S inner,
                      S outer,
                      Event e,
                      S& detector,
                      F& depth) {
    Point p1(0, 0), p2(0, 0);
    return check_for_hits(gen, model, inner, outer, e, detector, depth, p1, p2);
  }

  /// @param model acceptance model
  /// @param rx, ry coordinates of the emission point
  /// @param output parameter contains the lor of the event
  template <class RandomGenerator, class AcceptanceModel>
  short emit_event(RandomGenerator& gen,
                   AcceptanceModel& model,
                   F rx,
                   F ry,
                   F angle,
                   LOR& lor,
                   F& position) {

    typename Circle::Event e(rx, ry, angle);

    auto inner_secant = c_inner_.secant(e);
    auto outer_secant = c_outer_.secant(e);

    auto i_inner =
        c_inner_.section(c_inner_.angle(inner_secant.first), n_detectors_);
    auto i_outer =
        c_outer_.section(c_inner_.angle(outer_secant.first), n_detectors_);
    S detector1;
    F depth1;

    Point d1_p1(0, 0), d1_p2(0, 0);
    if (!check_for_hits(
            gen, model, i_inner, i_outer, e, detector1, depth1, d1_p1, d1_p2))
      return 0;

    i_inner =
        c_inner_.section(c_inner_.angle(inner_secant.second), n_detectors_);
    i_outer =
        c_outer_.section(c_inner_.angle(outer_secant.second), n_detectors_);
    S detector2;
    F depth2;
    Point d2_p1(0, 0), d2_p2(0, 0);
    if (!check_for_hits(
            gen, model, i_inner, i_outer, e, detector2, depth2, d2_p1, d2_p2))
      return 0;

    lor = LOR(detector1, detector2);

    if (lor.first == lor.second) {
      std::ostringstream msg;
      msg << __PRETTY_FUNCTION__ << " invalid LOR (" << lor.first << ", "
          << lor.second << ")";
      throw(msg.str());
    }

    Point origin(rx, ry);
    F length1 = origin.nearest_distance(d1_p1, d1_p2);
    length1 += depth1;
    F length2 = origin.nearest_distance(d2_p1, d2_p2);
    length2 += depth2;

    if (detector1 > detector2) {
      position = length1 - length2;
    } else {
      position = length2 - length1;
    }

#if DEBUG
    std::cerr << "position " << position << std::endl;
#endif
    return 2;
  }

  friend svg_ostream<F>& operator<<(svg_ostream<F>& svg, DetectorRing& dr) {
    svg << dr.c_outer_;
    svg << dr.c_inner_;

    for (auto detector = dr.begin(); detector != dr.end(); ++detector) {
      svg << *detector;
    }

    return svg;
  }

 private:
  Circle c_inner_;
  Circle c_outer_;
  S n_detectors_;
  S n_lors_;
  F fov_radius_;
  F radius_diff_;
};

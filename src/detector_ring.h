#pragma once

#include <map>

#include "random.h"
#include "detector.h"
#include "circle.h"
#include "svg_ostream.h"

/// Provides model for 2D ring of detectors
template <typename F = double>
class DetectorRing : public std::vector<Detector<F>> {
 public:
  typedef std::pair<size_t, size_t> LOR;
  typedef ::Circle<F> Circle;
  typedef ::Point<F> Point;
  typedef ::Detector<F> Detector;
  typedef ::Event<F> Event;

  /// @param n_detectors number of detectors on ring
  /// @param radius      radius of ring
  /// @param w_detector  width of single detector (along ring)
  /// @param h_detector  height/depth of single detector
  ////                   (perpendicular to ring)
  DetectorRing(size_t n_detectors, F radius, F w_detector, F h_detector)
      : c_inner_(radius),
        c_outer_(radius + h_detector),
        n_detectors_(n_detectors),
        n_lors_(n_detectors * (n_detectors + 1) / 2) {
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
  size_t lors() const { return n_lors_; }
  size_t detectors() const { return n_detectors_; }
  F fov_radius() const { return fov_radius_; }

  std::pair<int, int> pixel(F x, F y, F pixel_size) {
    F rx = x + fov_radius();
    F ry = y + fov_radius();
    return std::make_pair(static_cast<int>(floor(rx / pixel_size)),
                          static_cast<int>(floor(ry / pixel_size)));
  }

  template <class RandomGenerator, class AcceptanceModel>
  int check_for_hits(RandomGenerator& gen,
                     AcceptanceModel& model,
                     int inner,
                     int outer,
                     Event e,
                     int& detector) {

    // tells in which direction we got shorter modulo distance
    int step =
        ((n_detectors_ + inner - outer) % n_detectors_ >
         (n_detectors_ + outer - inner) % n_detectors_) ? 1 : n_detectors_ - 1;
    int end = (outer + step) % n_detectors_;
    for (int i = inner; i != end; i = (i + step) % n_detectors_) {
      auto points = (*this)[i].intersections(e);
      // check if we got 2 point intersection
      // then test the model against these points distance
      if (points.size() == 2) {
        if (model.deposition_depth(gen) < (points[1] - points[0]).length()) {
          detector = i;
          return 1;
        }
      }
    }
    return 0;
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
                   LOR& lor) {

    typename Circle::Event e(rx, ry, angle);

    auto inner_secant = c_inner_.secant(e);
    auto outer_secant = c_outer_.secant(e);

    auto i_inner =
        c_inner_.section(c_inner_.angle(inner_secant.first), n_detectors_);
    auto i_outer =
        c_outer_.section(c_inner_.angle(outer_secant.first), n_detectors_);
    int detector1;
    if (check_for_hits(gen, model, i_inner, i_outer, e, detector1) == 0)
      return 0;

    i_inner =
        c_inner_.section(c_inner_.angle(inner_secant.second), n_detectors_);
    i_outer =
        c_outer_.section(c_inner_.angle(outer_secant.second), n_detectors_);
    int detector2;
    if (check_for_hits(gen, model, i_inner, i_outer, e, detector2) == 0)
      return 0;

    lor.first = detector1;
    lor.second = detector2;
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
  size_t n_detectors_;
  size_t n_lors_;
  F fov_radius_;
};

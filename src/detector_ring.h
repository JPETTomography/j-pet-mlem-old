#pragma once

#include <map>

#include "random.h"
#include "detector.h"
#include "circle.h"
#include "svg_ostream.h"

/// Provides model for 2D ring of detectors
template <typename F = double>
class detector_ring : public std::vector<detector<F>> {
public:
  typedef std::pair<size_t, size_t> lor_type;
  typedef circle<F>   circle_type;
  typedef point<F>    point_type;
  typedef detector<F> detector_type;
  typedef event<F>    event_type;

  /// @param a_n_detectors number of detectors on ring
  /// @param radious       radious of ring
  /// @param w_detector    width of single detector (along ring)
  /// @param h_detector    height/depth of single detector
  ////                     (perpendicular to ring)
  detector_ring(size_t a_n_detectors, F radious, F w_detector, F h_detector)
  : c_inner(radious)
  , c_outer(radious+h_detector)
  , n_detectors(a_n_detectors)
  , n_lors( a_n_detectors * (a_n_detectors+1) / 2 )
  {
    if(radious    <= 0.)   throw("invalid radious");
    if(w_detector <= 0. ||
       h_detector <= 0.)   throw("invalid detector size");
    if(n_detectors % 4)    throw("number of detectors must be multiple of 4");

    fov_radius_ = radious / M_SQRT2;

    detector_type detector_base(h_detector, w_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base += point<>(radious + h_detector/2, 0.);

    // produce detector ring rotating base detector n times
    std::vector<detector<>> detector_ring;
    for (auto n = 0; n < n_detectors; ++n) {
      this->push_back( detector_base.rotated(2. * M_PI * n / n_detectors) );
    }
  }

  F radious()        const { return c_inner.radious(); }
  size_t lors()      const { return n_lors;            }
  size_t detectors() const { return n_detectors;       }
  F fov_radius()     const { return fov_radius_;       }

  std::pair<int,int> pixel(F x, F y, F pixel_size) {
    F rx = x + fov_radius();
    F ry = y + fov_radius();
    return std::make_pair(static_cast<int>( floor(rx / pixel_size) ),
                          static_cast<int>( floor(ry / pixel_size) ));
  }

  /// @param model acceptance model
  ///        (returns bool for call operator with given length)
  /// @param rx, ry coordinates of the emission point
  /// @param output parameter contains the lor of the event
  template <class RandomGenerator, class AcceptanceModel>
  short emit_event(RandomGenerator &gen,
                   AcceptanceModel &model,
                   F rx,
                   F ry,
                   F angle,
                   lor_type &lor) {

    typename circle_type::event_type e(rx, ry, angle);

    // secant for p and phi
    auto i_inner = c_inner.secant_sections(e, n_detectors);
    auto i_outer = c_outer.secant_sections(e, n_detectors);

    auto inner = i_inner.first;
    auto outer = i_outer.first;

    auto hits = 0;

#if COLLECT_INTERSECTIONS
    std::vector<point_type> ipoints;
#endif

    // process possible hit detectors on both sides
    for (auto side = 0; side < 2; ++side) {
      // starting from inner index
      // iterating to outer index
      auto i = inner;
      auto prev_i = i;

      // tells in which direction we got shorter modulo distance
      auto step = ((n_detectors+inner-outer) % n_detectors
                   >
                   (n_detectors+outer-inner) % n_detectors
                   ) ? 1 : -1;

#if SKIP_INTERSECTION
      (!(hits++) ? lor.first : lor.second) = i;
#else
      do {
        auto points = (*this)[i].intersections(e);
        // check if we got 2 point intersection
        // then test the model against these points distance

        if ( points.size() == 2 &&
             model( gen, (points[1]-points[0]).length() ) ) {

          hits++;
          (!side ? lor.first : lor.second) = i;
#if COLLECT_INTERSECTIONS
          for(auto &p: points) ipoints.push_back(p);
#endif
          break;
        }
        // step towards outer detector

        prev_i = i, i = (i + step + n_detectors) % n_detectors;

      } while (prev_i != outer); // loop over intersected detectors
#endif
      if (hits == 0) break;

      // switch side
      inner = i_inner.second;
      outer = i_outer.second;
    }

#if COLLECT_INTERSECTIONS
    if (hits >= 2)
      for(auto &p: ipoints) intersection_points.push_back(p);
#endif

    return hits;
  }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, detector_ring &dr) {
    svg << dr.c_outer;
    svg << dr.c_inner;

    for (auto detector = dr.begin(); detector != dr.end(); ++detector) {
      svg << *detector;
    }
#if COLLECT_INTERSECTIONS
    for (auto it = dr.intersection_points.begin(), p = *it; it != dr.intersection_points.end(); ++p, p = *it) {
      svg << "<circle cx=\"" << p.x << "\" cy=\"" << p.y << "\" r=\"0.002\"/>" << std::endl;
    }
#endif
    return svg;
  }

private:
  circle_type c_inner;
  circle_type c_outer;
  size_t n_detectors;
  size_t n_lors;
  F fov_radius_;
};

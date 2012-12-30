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
  /// @param radius       radius of ring
  /// @param w_detector    width of single detector (along ring)
  /// @param h_detector    height/depth of single detector
  ////                     (perpendicular to ring)
  detector_ring(size_t a_n_detectors, F radius, F w_detector, F h_detector)
  : c_inner(radius)
  , c_outer(radius+h_detector)
  , n_detectors(a_n_detectors)
  , n_lors( a_n_detectors * (a_n_detectors+1) / 2 ) {
    if(radius    <= 0.)   throw("invalid radius");
    if(w_detector <= 0. ||
       h_detector <= 0.)   throw("invalid detector size");
    if(n_detectors % 4)    throw("number of detectors must be multiple of 4");

    fov_radius_ = radius / M_SQRT2;

    detector_type detector_base(h_detector, w_detector);

    // move detector to the right edge of inner ring
    // along zero angle polar coordinate
    detector_base += point<>(radius + h_detector/2, 0.);

    // produce detector ring rotating base detector n times
    std::vector<detector<>> detector_ring;
    for (auto n = 0; n < n_detectors; ++n) {
      this->push_back( detector_base.rotated(2. * M_PI * n / n_detectors) );
    }
  }

  F radius()        const { return c_inner.radius(); }
  size_t lors()      const { return n_lors;            }
  size_t detectors() const { return n_detectors;       }
  F fov_radius()     const { return fov_radius_;       }

  std::pair<int,int> pixel(F x, F y, F pixel_size) {
    F rx = x + fov_radius();
    F ry = y + fov_radius();
    return std::make_pair(static_cast<int>( floor(rx / pixel_size) ),
                          static_cast<int>( floor(ry / pixel_size) ));
  }

  template <class RandomGenerator, class AcceptanceModel>
  int check_for_hits(RandomGenerator &gen,
                 AcceptanceModel &model,
                 int inner,
                 int outer,
                 event_type e, 
                 int &detector) {

    // tells in which direction we got shorter modulo distance
    int step = ((n_detectors+inner-outer) % n_detectors
                 >
                 (n_detectors+outer-inner) % n_detectors
                 ) ? 1 : n_detectors-1;
    int end=(outer+step)%n_detectors;
    int hits=0;
    int prev_i;
    int i=inner;
    for(int i=inner;i!=end;i=(i+step)%n_detectors) {
      auto points = (*this)[i].intersections(e);
      // check if we got 2 point intersection
      // then test the model against these points distance
      if ( points.size() == 2 ) {
        if(model.deposition_depth( gen) < (points[1]-points[0]).length() ) {
          detector=i;
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
  short emit_event(RandomGenerator &gen,
                   AcceptanceModel &model,
                   F rx,
                   F ry,
                   F angle,
                   lor_type &lor) {

    typename circle_type::event_type e(rx, ry, angle);

    auto inner_secant = c_inner.secant(e);
    auto outer_secant = c_outer.secant(e);


    auto i_inner = c_inner.section(c_inner.angle(inner_secant.first), n_detectors);
    auto i_outer = c_outer.section(c_inner.angle(outer_secant.first), n_detectors);
    int detector1;
    if(check_for_hits(gen,model,i_inner,i_outer,e,detector1)==0) return 0;

    i_inner = c_inner.section(c_inner.angle(inner_secant.second), n_detectors);
    i_outer = c_outer.section(c_inner.angle(outer_secant.second), n_detectors);
    int detector2;
    if(check_for_hits(gen,model,i_inner,i_outer,e,detector2)==0) return 0;

    lor.first=detector1;
    lor.second=detector2;
    return 2;
  }

  friend svg_ostream<F> & operator << (svg_ostream<F> &svg, detector_ring &dr) {
    svg << dr.c_outer;
    svg << dr.c_inner;

    for (auto detector = dr.begin(); detector != dr.end(); ++detector) {
      svg << *detector;
    }

    return svg;
  }

private:
  circle_type c_inner;
  circle_type c_outer;
  size_t n_detectors;
  size_t n_lors;
  F fov_radius_;
};

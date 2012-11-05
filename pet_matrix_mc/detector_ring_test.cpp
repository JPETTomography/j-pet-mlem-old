#include <iostream>
#include <fstream>

#include <catch.hpp>

#include "model.h"
#include "detector_ring.h"

template<typename T, typename I>
bool is_in(T val, I it, int n) {

  for(int i=0;i<n;++i) {
    if(it[i]==val) return true;
  }

  return false;
}

TEST_CASE("detector_ring/math", "detector ring test") {

  std::ifstream  in("detector_ring.test");
  if(in) {
    double r,w,h;
    int n_detectors;
    int n_events;
    int n_pixels=128;

    in>>r>>n_detectors>>w>>h>>n_events;;

    double s_pixel = r / (n_pixels*sqrt(2.0));

    detector_ring<> ring(n_detectors,n_pixels,s_pixel,r,w,h);

    for(int i_event=0;i_event<n_events;++i_event) {
      double x,y,phi;
      in>>x>>y>>phi;

      decltype(ring)::event_type event(x,y,phi);

      in>>n_detectors;

      int detector[n_detectors];
      for(int i=0;i<n_detectors;i++) {
        in>>detector[i];
        detector[i]--; //mathematica counts positions from 1
      }

      for(int i=0;i<n_detectors;i++) {
        double x,y;
        in>>x>>y;
        decltype(ring)::point_type p1(x,y);
        in>>x>>y;
        decltype(ring)::point_type p2(x,y);

        auto inters = ring[detector[i]].intersections(event);
        CHECK(inters.size()==2);

        double tol=1e-14;

        bool first_to_first =compare(p1,inters[0],tol) &&
          compare(p2,inters[1],tol);
        bool first_to_second=compare(p1,inters[1],tol) &&
          compare(p2,inters[0],tol);

        CHECK( (first_to_first || first_to_second)== true);
      }

      //this is not yet a complete tests....
      decltype(ring)::lor_type lor;
      always_accept<> model;
      auto hits=ring.emit_event(model,x,y,phi,lor);

      if(hits>=2) {
        CHECK(is_in(lor.first ,detector,n_detectors)==true);
        CHECK(is_in(lor.second,detector,n_detectors)==true);
      }
    }
  }
  else {
    WARN("cannot open file `detector_ring.test'");
  }
}

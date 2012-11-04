#include<iostream>
#include<fstream>

#include<catch.hpp>

#include"detector_ring.h"

TEST_CASE("detector_ring/math", "detector ring test") {

  std::ifstream  in("detector_ring.test");
  double r,w,h;
  int n_detectors;
  int n_events;
  int n_pixels=128;

  in>>r>>n_detectors>>w>>h>>n_events;;

  //  std::cerr<<r<<" "<<n_detectors<<" "<<w<<" "<<h<<" "<<n_events<<std::endl;
  double s_pixel=r/(n_pixels*sqrt(2.0));

  detector_ring<> ring(n_detectors,n_pixels,s_pixel,r,w,h);

  for(int i_event=0;i_event<n_events;++i_event) {
    double x,y,phi;
    in>>x>>y>>phi;
    
    decltype(ring)::event_type event(x,y,phi);

    in>>n_detectors;
    // std::cerr<<n_detectors<<std::endl;
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
 
  }
}

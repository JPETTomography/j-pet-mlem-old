#include<iostream>
#include<fstream>

#include <catch.hpp>

#include "polygon.h"

TEST_CASE("polygon/intersection", "polygon intersection") {
  polygon<> p;
  p.push_back( {1., 1.} );
  p.push_back( {2., 1.} );
  p.push_back( {2., 2.} );
  p.push_back( {1., 2.} );

  decltype(p)::event_type e1(0.,  .5, M_PI_4);
  decltype(p)::event_type e2(0., 1.5, M_PI_4);
  decltype(p)::event_type e3(0., 1.5, 0.);

  CHECK( true  == p.intersects(e1) );
  CHECK( false == p.intersects(e2) );
  CHECK( true  == p.intersects(e3) );

  SECTION("polygon/intersection/points", "intersection points") {
    auto i1 = p.intersections(e1);

    REQUIRE( i1.size() == 2 );
    CHECK( std::min(i1[0].x, i1[1].x) == Approx(1.)  );
    CHECK( std::max(i1[0].x, i1[1].x) == Approx(1.5) );

    CHECK( std::min(i1[0].y, i1[1].y) == Approx(1.5) );
    CHECK( std::max(i1[0].y, i1[1].y) == Approx(2.)  );

    auto i3 = p.intersections(e3);

    REQUIRE( i3.size() == 2 );
    CHECK( std::min(i3[0].x, i3[1].x) == Approx(1.)  );
    CHECK( std::max(i3[0].x, i3[1].x) == Approx(2.)  );

    CHECK( std::min(i3[0].y, i3[1].y) == Approx(1.5) );
    CHECK( std::max(i3[0].y, i3[1].y) == Approx(1.5) );
  }
}

TEST_CASE("polygon/intersection/math","rectangle intersections from mathematica") {
  std::ifstream in("polygon.test");

  if(!in) {
    WARN("cannot open file `polygon.test'");
    return;
  }

  int n_events;
  in>>n_events;
  std::cerr<<n_events<<std::endl;

  polygon<> poly;

  for(int i=0;i<4;++i) {
    double x,y;
    in>>x>>y;
    decltype(poly)::point_type p(x,y);
    poly.push_back(p);
  }

  for(int i=0;i<n_events;i++) {
    double x,y,phi;
    in>>x>>y>>phi;

    double a,b,c;
    in>>a>>b>>c;

    int n_itersections;
    in >> n_itersections;

    decltype(poly)::event_type event(x,y,phi);
    bool intersects=n_itersections>0;
    const double tol=1e-14;
    CHECKED_IF(poly.intersects(event)==intersects) {
      auto intersections=poly.intersections(event);
      CHECKED_IF(intersections.size()==n_itersections) {

        if(n_itersections>0) {
          double ix,iy;
          decltype(poly)::intersections_type m_intersections;
          for(int j=0;j<n_itersections;++j) {
            in>>ix>>iy;
            decltype(poly)::point_type p(ix,iy);
            m_intersections.push_back(p);
          }

          if(n_itersections==1) {
            CHECK(compare(intersections[0],m_intersections[0],tol)==true);
          } else {
            bool first_to_first =
              compare(intersections[0],m_intersections[0],tol) &&
              compare(intersections[1],m_intersections[1],tol);

            bool first_to_second =
              compare(intersections[0],m_intersections[1],tol) &&
              compare(intersections[1],m_intersections[0],tol);

            CHECK( (first_to_first || first_to_second)==true);
          }
        }
      }
    }
  }
}

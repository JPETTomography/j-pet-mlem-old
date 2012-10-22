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

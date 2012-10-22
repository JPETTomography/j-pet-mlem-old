#include <catch.hpp>

#include "circle.h"

TEST_CASE("circle/init", "circle initialization") {
  circle<> c1(1.);

  CHECK( c1.radious()  == 1. );
  CHECK( c1.radious2() == 1. );

  circle<> c2(std::sqrt(2.));

  CHECK( c2.radious()  == std::sqrt(2.) ); // exact!
  CHECK( c2.radious2() == Approx(2.)    );
}

TEST_CASE("circle/secant", "circle secant") {
  circle<> c(1);

  SECTION("angle-0", "0 degrees from (0, 0)") {
    decltype(c)::event_type zero(0., 0., 0.);
    auto s = c.secant(zero);

    CHECK( std::min(s.first.x, s.second.x) == Approx(-1.) );
    CHECK( std::max(s.first.x, s.second.x) == Approx( 1.) );

    CHECK( s.first.y  ==  0. );
    CHECK( s.second.y ==  0. );

    auto a = c.secant_angles(zero);
    if (a.first  == Approx(-M_PI)) a.first  += 2. * M_PI;
    if (a.second == Approx(-M_PI)) a.second += 2. * M_PI;

    CHECK( std::min(a.first, a.second) == Approx(0.) );
    CHECK( std::max(a.first, a.second) == Approx(M_PI) );
  }
  SECTION("angle-90", "90 degrees from (0, 0)") {
    decltype(c)::event_type zero90(0., 0., M_PI_2);
    auto s = c.secant(zero90);

    CHECK( s.first.x  == Approx(0.) );
    CHECK( s.second.x == Approx(0.) );

    CHECK( std::min(s.first.y, s.second.y) == -1. );
    CHECK( std::max(s.first.y, s.second.y) ==  1. );

    auto a = c.secant_angles(zero90);
    if (a.first  == Approx(-M_PI)) a.first  += M_2_PI;
    if (a.second == Approx(-M_PI)) a.second += M_2_PI;

    CHECK( std::min(a.first, a.second) == Approx(-M_PI_2) );
    CHECK( std::max(a.first, a.second) == Approx( M_PI_2) );
  }
#if 0
  SECTION("angle-45", "45 degrees from (1, 0)") {
    decltype(c)::event_type xone45(1., 0., M_PI_4);
    auto s = c.secant(xone45);

    CHECK( std::min(s.first.x, s.second.x) == Approx( 0.) );
    CHECK( std::max(s.first.x, s.second.x) == xone45.x    );

    CHECK( std::min(s.first.y, s.second.y) == Approx(-1.) );
    CHECK( std::max(s.first.y, s.second.y) == xone45.y    );
  }
#endif
}

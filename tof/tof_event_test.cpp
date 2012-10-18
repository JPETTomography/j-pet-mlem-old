#include <catch.hpp>

#include <cmath>

#include "tof_event.h"

TEST_CASE("tof_event", "TOF event") {

  ToF_Event_2D<double> event(200.0, 150.0, 1.0);

  REQUIRE( 200 == event.z()   );
  REQUIRE( 150.0 == event.y() );
  REQUIRE( 1.0 == event.tan() );
}

#if 0
TEST_CASE("tof_event/from_ps", "TOF event from PS") {

  ToF_Event_2D<double> event;
  ToF_Event_2D<double>::fromPS(event, 350.0, -350.0, 0.0, 350.0);

  REQUIRE( 0.0 == event.z()   );
  REQUIRE( 0.0 == event.y()   );
  REQUIRE( 1.0 == event.tan() );
}
#endif

#include <catch.hpp>

#include <cmath>

#include "tof_event.h"
#include "tof_detector.h"
#include "phantom.h"
#include "tof_monte_carlo.h"

typedef ToF_Detector_2D<double> detector_t;
typedef ToF_Event_2D<double> event_t;

TEST_CASE("monte_carlo", "Monte-Carlo") {

  detector_t detector(350, 500);
  detector.set_sigma(11, 32);

  ToF_Monte_Carlo<double, detector_t> mc(detector);
  mc.seed(544445);

  SECTION("add_noise_to_track", "") {

    ToF_Track_2D<double> track_in(-200, 23, 117);
    ToF_Track_2D<double> track_out = mc.add_noise(track_in);

    REQUIRE(track_in.z_up() != track_out.z_up());
    REQUIRE(track_in.z_dn() != track_out.z_dn());
    REQUIRE(track_in.dl() != track_out.dl());
  }

  SECTION("add_noise_to_event", "") {

    ToF_Event_2D<double> event_in(-10, 23, -.75);
    ToF_Event_2D<double> event_out = mc.add_noise(event_in);

    REQUIRE(event_in.z() != event_out.z());
    REQUIRE(event_in.y() != event_out.y());
    REQUIRE(event_in.tan() != event_out.tan());
  }
}

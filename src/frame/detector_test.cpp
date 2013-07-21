#include "catch.hpp"

typedef float FLOAT;

#include"detector.h"

TEST_CASE("EventImage/Angle", "Angle") {

  EventImageAngle<FLOAT> event(20.0, 100.0, 18 * Degree);

  REQUIRE(event.y() == Approx(20.0));
  REQUIRE(event.z() == Approx(100.0));
  REQUIRE(event.theta() == Approx(18 * Degree));
}

TEST_CASE("EventImage/Tan", "Tangent") {

  EventImageTan<FLOAT> event(20.0, 100.0, 0.75);

  REQUIRE(event.y() == Approx(20.0));
  REQUIRE(event.z() == Approx(100.0));
  REQUIRE(event.tan() == Approx(0.75));
}

TEST_CASE("EventImage/Convert", "Conversion") {

  EventImageTan<FLOAT> eventTan(20.0, 100.0, ::tan(18 * Degree));
  EventImageAngle<FLOAT> eventAngle = EventImageTanToAngle(eventTan);

  REQUIRE(eventAngle.y() == Approx(20.0));
  REQUIRE(eventAngle.z() == Approx(100.0));
  REQUIRE(eventAngle.theta() == Approx(18 * Degree));

  EventImageTan<FLOAT> eventTan2 = EventImageAngleToTan(eventAngle);

  REQUIRE(eventTan2.y() == Approx(20.0));
  REQUIRE(eventTan2.z() == Approx(100.0));
  REQUIRE(eventTan2.tan() == Approx(tan(18 * Degree)));
}

TEST_CASE("EventDetector/Construct", "Construc") {

  EventDetector<FLOAT> event(-20.0, 100.0, 400);

  REQUIRE(event.z_up() == Approx(-20.0));
  REQUIRE(event.z_dn() == Approx(100.0));
  REQUIRE(event.dl() == Approx(400.0));
}

TEST_CASE("Detector/Construct", "Construct") {

  Detector<FLOAT> detector(450.0, 300.0);
  REQUIRE(detector.R() == Approx(450.0));
  REQUIRE(detector.L() == Approx(300.0));
}

TEST_CASE("Detector/Conversions/ToEventImage", "Convert") {

  Detector<FLOAT> detector(450.0, 300.0);

  EventDetector<FLOAT> eventDetector(-70.0, 50.0, 100.0);

  EventImageTan<FLOAT> eventImageTan =
      detector.EventDetectorToImageTan(eventDetector);

  REQUIRE(eventImageTan.y() == Approx(-49.5613950341317));
  REQUIRE(eventImageTan.z() == Approx(-3.39181399544910));
  REQUIRE(eventImageTan.tan() == Approx(-0.13333333333333));

  EventImageAngle<FLOAT> eventImageAngle =
      detector.EventDetectorToImageAngle(eventDetector);

  REQUIRE(eventImageAngle.y() == Approx(-49.5613950341317));
  REQUIRE(eventImageAngle.z() == Approx(-3.39181399544910));
  REQUIRE(eventImageAngle.theta() == Approx(-0.132551532296674));
}

TEST_CASE("Detector/Conversions/ToEventDetector", "Convert") {
  Detector<FLOAT> detector(450.0, 300.0);

  EventImageAngle<FLOAT> eventImageAngle(100, 25, -18.0 * Degree);
  EventImageTan<FLOAT> eventImageTan = EventImageAngleToTan(eventImageAngle);

  EventDetector<FLOAT> eventDetectorTan =
      detector.EventImageTanToDetector(eventImageTan);
  REQUIRE(eventDetectorTan.z_up() == Approx(-88.72189368151721));
  REQUIRE(eventDetectorTan.z_dn() == Approx(203.7058329280985));
  REQUIRE(eventDetectorTan.dl() == Approx(-210.2924448476534));

  EventDetector<FLOAT> eventDetectorAngle =
      detector.EventImageAngleToDetector(eventImageAngle);

  REQUIRE(eventDetectorAngle.z_up() == Approx(-88.72189368151721));
  REQUIRE(eventDetectorAngle.z_dn() == Approx(203.7058329280985));
  REQUIRE(eventDetectorAngle.dl() == Approx(-210.2924448476534));
}

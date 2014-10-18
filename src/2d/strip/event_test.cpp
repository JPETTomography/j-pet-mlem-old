#include "catch.hpp"
#include <cmath>

#include "strip_detector.h"

const double degree = M_PI / 180.0;

TEST_CASE("strip/event/conversions1") {
  StripDetector<> detector(450.0, 200.0, 200, 200, 5.0, 5.0, 10, 63);

  ImageSpaceEventAngle<> img_angle(10.0, 20.0, 7.0 * degree);
  ImageSpaceEventTan<> img_tan = img_angle.to_tan();
  CHECK(img_tan.y == Approx(10.0).epsilon(1e-13));
  CHECK(img_tan.z == Approx(20.0).epsilon(1e-13));
  CHECK(img_tan.tan == Approx(0.1227845609029046).epsilon(1e-13));

  Event<> proj = detector.to_projection_space_tan(img_tan);

  CHECK(proj.z_u == Approx(74.02520679727803).epsilon(1e-13));
  CHECK(proj.z_d == Approx(-36.480898015336116).epsilon(1e-13));
  CHECK(proj.dl == Approx(-20.15019650917697).epsilon(1e-13));

  ImageSpaceEventAngle<> re_img_angle =
      detector.from_projection_space_angle(proj);

  CHECK(re_img_angle.y == Approx(img_angle.y).epsilon(1e-13));
  CHECK(re_img_angle.z == Approx(img_angle.z).epsilon(1e-13));
  CHECK(re_img_angle.angle == Approx(img_angle.angle).epsilon(1e-13));
}

TEST_CASE("strip/event/conversions2") {
  StripDetector<> detector(450.0, 200.0, 200, 200, 5.0, 5.0, 10, 63);

  ImageSpaceEventAngle<> img_angle(-10.0, 37.0, -5.0 * degree);
  ImageSpaceEventTan<> img_tan = img_angle.to_tan();
  CHECK(img_tan.y == Approx(-10.0).epsilon(1e-13));
  CHECK(img_tan.z == Approx(37.0).epsilon(1e-13));
  CHECK(img_tan.tan == Approx(-0.08748866352592401).epsilon(1e-13));

  Event<> proj = detector.to_projection_space_tan(img_tan);

  CHECK(proj.z_u == Approx(-3.244785221925042).epsilon(1e-13));
  CHECK(proj.z_d == Approx(75.49501195140655).epsilon(1e-13));
  CHECK(proj.dl == Approx(20.076396750866948).epsilon(1e-13));

  ImageSpaceEventAngle<> re_img_angle =
      detector.from_projection_space_angle(proj);

  CHECK(re_img_angle.y == Approx(img_angle.y).epsilon(1e-13));
  CHECK(re_img_angle.z == Approx(img_angle.z).epsilon(1e-13));
  CHECK(re_img_angle.angle == Approx(img_angle.angle).epsilon(1e-13));
}

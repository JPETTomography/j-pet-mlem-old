#include "util/test.h"
#include <cmath>

#include "detector.h"

using namespace PET2D;
using namespace PET2D::Strip;

const double degree = M_PI / 180.0;

TEST("2d/strip/event/conversions1") {
  Detector<double, short> detector(450.0, 200.0, 200, 200, 5.0, 5.0, 10, 63);

  ImageSpaceEventAngle<double> img_angle(10.0, 20.0, 7.0 * degree);
  ImageSpaceEventTan<double> img_tan = img_angle.to_tan();
  CHECK(img_tan.y == 10._e13);
  CHECK(img_tan.z == 20._e13);
  CHECK(img_tan.tan == 0.1227845609029046_e13);

  Event<double> proj = detector.to_projection_space_tan(img_tan);

  CHECK(proj.z_u == 74.02520679727803_e13);
  CHECK(proj.z_d == -36.480898015336116_e13);
  CHECK(proj.dl == -20.15019650917697_e13);

  ImageSpaceEventAngle<double> re_img_angle =
      detector.from_projection_space_angle(proj);

  CHECK(re_img_angle.y == Approx(img_angle.y).epsilon(1e-13));
  CHECK(re_img_angle.z == Approx(img_angle.z).epsilon(1e-13));
  CHECK(re_img_angle.angle == Approx(img_angle.angle).epsilon(1e-13));
}

TEST("2d/strip/event/conversions2") {
  Detector<double, short> detector(450.0, 200.0, 200, 200, 5.0, 5.0, 10, 63);

  ImageSpaceEventAngle<double> img_angle(-10.0, 37.0, -5.0 * degree);
  ImageSpaceEventTan<double> img_tan = img_angle.to_tan();
  CHECK(img_tan.y == -10._e13);
  CHECK(img_tan.z == 37._e13);
  CHECK(img_tan.tan == -0.08748866352592401_e13);

  Event<double> proj = detector.to_projection_space_tan(img_tan);

  CHECK(proj.z_u == -3.244785221925042_e13);
  CHECK(proj.z_d == 75.49501195140655_e13);
  CHECK(proj.dl == 20.076396750866948_e13);

  ImageSpaceEventAngle<double> re_img_angle =
      detector.from_projection_space_angle(proj);

  CHECK(re_img_angle.y == Approx(img_angle.y).epsilon(1e-13));
  CHECK(re_img_angle.z == Approx(img_angle.z).epsilon(1e-13));
  CHECK(re_img_angle.angle == Approx(img_angle.angle).epsilon(1e-13));
}

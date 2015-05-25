#include "util/test.h"

#include "frame.h"

using Point = PET3D::Point<float>;
using Point2D = PET2D::Point<float>;

TEST("PET3D/geometry/frame" ) {

  Frame<float> frame(Point(-0.205, 0.0, 0.0), Point(0.205, 0.0, 0.0));

  Point p(1, 2, 3);

  auto proj = frame.project(p);
  REQUIRE(proj.x == 3.0_e7);
  REQUIRE(proj.y == -1.0_e7);
  REQUIRE(frame.distance(p)==2.0_e7);


}


TEST("PET3D/geometry/frame/45degree" ) {

  Frame<float> frame(Point(-0.205, -0.205, 0.0), Point(0.205, 0.205, 0.0));

  Point p(1, 2, 3);

  auto proj = frame.project(p);
  REQUIRE(proj.x == 3.0_e7);
  REQUIRE(proj.y == -2.12132034_e7);
  REQUIRE(frame.distance(p)==0.707106781187_e7);

}

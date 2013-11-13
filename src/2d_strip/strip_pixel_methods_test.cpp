#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "strip_detector.h"

typedef StripDetector<double>::Pixel Pixel;
typedef StripDetector<double>::Point Point;

TEST_CASE("strip/pixel/location") {

  StripDetector<double> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  Pixel p = detector.pixel_location(0.0, 0.0);

  CHECK(p.first == 100);
  CHECK(p.second == 100);

  // test boundary points
  p = detector.pixel_location(500.0, -500.0);

  CHECK(p.first == 0);
  CHECK(p.second == 0);

  p = detector.pixel_location(500.0, 500.0);

  CHECK(p.first == 0);
  CHECK(p.second == 200);

  p = detector.pixel_location(-500.0, -500.0);

  CHECK(p.first == 200);
  CHECK(p.second == 0);

  p = detector.pixel_location(-500.0, 500.0);

  CHECK(p.first == 200);
  CHECK(p.second == 200);
}

TEST_CASE("strip/pixel/center") {

  // space->image_space  y: [R,-R] ->[0,n_pixels_y], z:[-L/2,L/2] ->
  // [0,n_pixels_z]
  StripDetector<double> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  // test middle point pixel center
  Pixel p = detector.pixel_location(0.0, 0.0);
  Point pc = detector.pixel_center(p.first, p.second);

  CHECK(pc.first == -2.5);
  CHECK(pc.second == 2.5);

  // test -y,+z
  p = detector.pixel_location(-6.0, 3.0);
  pc = detector.pixel_center(p.first, p.second);

  CHECK(pc.first == -7.5);
  CHECK(pc.second == 2.5);

  // test +y,+z
  p = detector.pixel_location(6.0, 3.0);
  pc = detector.pixel_center(p.first, p.second);

  CHECK(pc.first == 7.5);
  CHECK(pc.second == 2.5);

  // test +y,-z
  p = detector.pixel_location(6.0, -3.0);
  pc = detector.pixel_center(p.first, p.second);

  CHECK(pc.first == 7.5);
  CHECK(pc.second == -2.5);

  // test -y,-z
  p = detector.pixel_location(-6.0, -3.0);
  pc = detector.pixel_center(p.first, p.second);

  CHECK(pc.first == -7.5);
  CHECK(pc.second == -2.5);
}

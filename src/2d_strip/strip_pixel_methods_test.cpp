#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include"reconstruction.h"

typedef Reconstruction<double>::Pixel Pixel;
typedef Reconstruction<double>::Point Point;

TEST_CASE("Strip pixel location",
          "Reconstruction<>::pixel_location method test") {

  Reconstruction<double> reconstructor(1, 500, 1000, 200, 5, 10, 63);

  Pixel p = reconstructor.pixel_location(0.0, 0.0);

  CHECK(p.first == 100);
  CHECK(p.second == 100);

  // test boundary points
  p = reconstructor.pixel_location(500.0, -500.0);

  CHECK(p.first == 0);
  CHECK(p.second == 0);

  p = reconstructor.pixel_location(500.0, 500.0);

  CHECK(p.first == 0);
  CHECK(p.second == 200);

  p = reconstructor.pixel_location(-500.0, -500.0);

  CHECK(p.first == 200);
  CHECK(p.second == 0);

  p = reconstructor.pixel_location(-500.0, 500.0);

  CHECK(p.first == 200);
  CHECK(p.second == 200);
}

TEST_CASE("Strip pixel center", "Reconstruction<>::pixel_center method test") {

  // space->image_space  y: [R,-R] ->[0,n_pixels_y], z:[-L/2,L/2] ->
  // [0,n_pixels_z]
  Reconstruction<double> reconstructor(1, 500, 1000, 200, 5, 10, 63);

  // test middle point pixel center
  Pixel p = reconstructor.pixel_location(0.0, 0.0);
  Point pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == -2.5);
  CHECK(pc.second == 2.5);

  // test -y,+z
  p = reconstructor.pixel_location(-6.0, 3.0);
  pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == -7.5);
  CHECK(pc.second == 2.5);

  // test +y,+z
  p = reconstructor.pixel_location(6.0, 3.0);
  pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == 7.5);
  CHECK(pc.second == 2.5);

  // test +y,-z
  p = reconstructor.pixel_location(6.0, -3.0);
  pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == 7.5);
  CHECK(pc.second == -2.5);

  // test -y,-z
  p = reconstructor.pixel_location(-6.0, -3.0);
  pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == -7.5);
  CHECK(pc.second == -2.5);
}

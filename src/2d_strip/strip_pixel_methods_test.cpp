#include "catch.hpp"
#include <cmath>
#include <utility>

#include"reconstruction.h"

typedef Reconstruction<double>::Pixel Pixel;
typedef Reconstruction<double>::Point Point;



TEST_CASE("strip_pixel_locations", "strip_pixel_methods_test") {
Reconstruction<double> reconstructor(1, 500, 1000, 200, 5, 10, 63);

  // test middle point
  Pixel p = reconstructor.pixel_location(0.0f, 0.0f);

  CHECK(p.first  == 100);
  CHECK(p.second == 100);

  // test boundary points
  p = reconstructor.pixel_location(500.0f, -500.0f);

  CHECK(p.first == 0);
  CHECK(p.second == 0);

  p = reconstructor.pixel_location(500.0f, 500.0f);

  CHECK(p.first == 0);
  CHECK(p.second == 200);

  p = reconstructor.pixel_location(-500.0f, -500.0f);

  CHECK(p.first == 200);
  CHECK(p.second == 0);

  p = reconstructor.pixel_location(-500.0f, 500.0f);

  CHECK(p.first == 200);
  CHECK(p.second == 200);

}



TEST_CASE("strip_pixel_centers", "strip_pixel_center_methods_test") {
Reconstruction<double> reconstructor(1, 500, 1000, 200, 5, 10, 63);
  // test middle point pixel center

  Pixel p =  reconstructor.pixel_location(0.0f, 0.0f);
  Point pc = reconstructor.pixel_center(p.first, p.second);

  CHECK(pc.first == 2.5);
  CHECK(pc.second == 2.5);


}


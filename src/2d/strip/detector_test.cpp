#include "util/test.h"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "detector.h"

using Detector = PET2D::Strip::Detector<>;
using Point = PET2D::Point<>;
using Pixel = PET2D::Pixel<>;

TEST("2d/strip/detector/pixel_at") {

  Detector d(500, 1000, 200, 200, 5, 5, 10, 63);

  SECTION("detector_center") {
    CHECK(d.pixel_at({ 0., 0. }) == Pixel(100, 100));
    CHECK(d.pixel_at({ -0.00001, 0.00001 }) == Pixel(99, 99));
  }

  SECTION("boundaries") {
    CHECK(d.pixel_at({ -500., 500. }) == Pixel(0, 0));      // TL
    CHECK(d.pixel_at({ 500., 500. }) == Pixel(200, 0));     // TR
    CHECK(d.pixel_at({ -500., -500. }) == Pixel(0, 200));   // BL
    CHECK(d.pixel_at({ 500., -500. }) == Pixel(200, 200));  // BR
  }

  SECTION("close_to_boundaries") {
    CHECK(d.pixel_at({ -499.99999, 499.99999 }) == Pixel(0, 0));      // TL
    CHECK(d.pixel_at({ 499.99999, 499.99999 }) == Pixel(199, 0));     // TR
    CHECK(d.pixel_at({ -499.99999, -499.99999 }) == Pixel(0, 199));   // BL
    CHECK(d.pixel_at({ 499.99999, -499.99999 }) == Pixel(199, 199));  // BR
  }
}

TEST("2d/strip/detector/pixel_center") {

  // space->image_space  y: [  R,  -R ] -> [0, n_pixels_y],
  //                     z: [-L/2, L/2] -> [0, n_pixels_z]

  SECTION("pixel_locations") {
    Detector d(100, 200, 200, 200, 10, 63);
    CHECK(d.pixel_center({ 100, 100 }) == Point(0.5, -0.5));
    CHECK(d.pixel_center({ 99, 99 }) == Point(-0.5, 0.5));
    CHECK(d.pixel_center({ 0, 0 }) == Point(-99.5, 99.5));
  }

  SECTION("reversibility_1") {
    Detector d(100, 200, 200, 200, 10, 63);
    CHECK(d.pixel_at(d.pixel_center({ 100, 100 })) == Pixel(100, 100));
    CHECK(d.pixel_at(d.pixel_center({ 99, 99 })) == Pixel(99, 99));
    CHECK(d.pixel_at(d.pixel_center({ 0, 0 })) == Pixel(0, 0));
  }

  SECTION("reversibility_2") {
    Detector d(500, 1000, 200, 200, 5, 5, 10, 63);  // middle
    CHECK(d.pixel_center(d.pixel_at({ 0., 0. })) == Point(2.5, -2.5));
    CHECK(d.pixel_center(d.pixel_at({ -6., 3. })) == Point(-7.5, 2.5));    // TL
    CHECK(d.pixel_center(d.pixel_at({ 6., 3. })) == Point(7.5, 2.5));      // TR
    CHECK(d.pixel_center(d.pixel_at({ -6., -3. })) == Point(-7.5, -2.5));  // BL
    CHECK(d.pixel_center(d.pixel_at({ 6., -3. })) == Point(7.5, -2.5));    // BR
  }
}

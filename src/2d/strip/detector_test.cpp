#include "util/test.h"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "detector.h"

using namespace PET2D;
using namespace PET2D::Strip;

TEST("strip/detector/pixel_at") {

  Detector<> d(500, 1000, 200, 200, 5, 5, 10, 63);

  CHECK(d.pixel_at({ 0.0, 0.0 }) == Pixel<>(100, 100));

  // test boundary points
  CHECK(d.pixel_at({ 500.0, -500.0 }) == Pixel<>(0, 0));
  CHECK(d.pixel_at({ 500.0, 500.0 }) == Pixel<>(0, 200));
  CHECK(d.pixel_at({ -500.0, -500.0 }) == Pixel<>(200, 0));
  CHECK(d.pixel_at({ -500.0, 500.0 }) == Pixel<>(200, 200));
}

TEST("strip/detector/pixel_center") {

  // space->image_space  y: [  R,  -R ] -> [0, n_pixels_y],
  //                     z: [-L/2, L/2] -> [0, n_pixels_z]
  Detector<> d(500, 1000, 200, 200, 5, 5, 10, 63);

  // test middle point pixel center
  CHECK(d.pixel_center(d.pixel_at({ 0.0, 0.0 })) == Point<>(-2.5, 2.5));
  // test -y,+z
  CHECK(d.pixel_center(d.pixel_at({ -6.0, 3.0 })) == Point<>(-7.5, 2.5));
  // test +y,+z
  CHECK(d.pixel_center(d.pixel_at({ 6.0, 3.0 })) == Point<>(7.5, 2.5));
  // test +y,-z
  CHECK(d.pixel_center(d.pixel_at({ 6.0, -3.0 })) == Point<>(7.5, -2.5));
  // test -y,-z
  CHECK(d.pixel_center(d.pixel_at({ -6.0, -3.0 })) == Point<>(-7.5, -2.5));
}

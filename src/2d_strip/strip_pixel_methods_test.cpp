#include "catch.hpp"
#include <cmath>
#include <utility>

typedef std::pair<int, int> Pixel;
typedef std::pair<float, float> Point;

template <typename T> Pixel pixel_location(T y, T z) {
  T R_distance = 500.0f;
  T pixel_size = 5.0f;
  return Pixel(std::floor((R_distance - y) / pixel_size),
               std::floor((R_distance + z) / pixel_size));
}

template <typename T> Point pixel_center(T y, T z) {
  T R_distance = 500.0f;
  T pixel_size = 5.0f;
  return Point(
      (std::floor(R_distance - (y) * pixel_size)) + (T(0.5) * pixel_size),
      (std::floor((z) * pixel_size - R_distance)) + (T(0.5) * pixel_size));
}

TEST_CASE("strip_pixel_locations", "strip_pixel_methods_test") {

  // test middle point
  Pixel p = pixel_location<float>(0.0f, 0.0f);

  CHECK(p.first == 100);
  CHECK(p.second == 100);

  // test boundary points
  p = pixel_location<float>(500.0f, -500.0f);

  CHECK(p.first == 0);
  CHECK(p.second == 0);

  p = pixel_location<float>(500.0f, 500.0f);

  CHECK(p.first == 0);
  CHECK(p.second == 200);

  p = pixel_location<float>(-500.0f, -500.0f);

  CHECK(p.first == 200);
  CHECK(p.second == 0);

  p = pixel_location<float>(-500.0f, 500.0f);

  CHECK(p.first == 200);
  CHECK(p.second == 200);

  // CHECK(vec_o[1] == Approx(-2100.0));
  // CHECK(vec_o[2] == Approx(-1800 * std::sqrt(2)));
}

TEST_CASE("strip_pixel_centers", "strip_pixel_center_methods_test") {

  // test middle point pixel center
  Pixel p = pixel_location<float>(0.0f, 0.0f);
  Point pc = pixel_center<float>(p.first, p.second);

  CHECK(pc.first == 2.5);
  CHECK(pc.second == 2.5);

  // CHECK(vec_o[1] == Approx(-2100.0));
  // CHECK(vec_o[2] == Approx(-1800 * std::sqrt(2)));
}

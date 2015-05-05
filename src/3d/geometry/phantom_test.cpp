#include <random>

#include "util/test.h"

#include "phantom.h"

TEST("PET3D/Geometry/CylinderRegion") {
  using F = float;
  using Region =
      CylinderRegion<F, std::mt19937, PET3D::spherical_distribution<F>>;
  using Point = Region::Point;

  Region region(
      2.0f, 3.0f, 1.0f, PET3D::spherical_distribution<F>(-M_PI / 3, M_PI / 3));

  REQUIRE(region.volume() == Approx(4.0 * M_PI * 3.0).epsilon(1e-7));
  REQUIRE(region.intensity == 1.0_e7);

  Point p1(1.2f, 0.1f, 1.4f);
  REQUIRE(region.in(p1));
  Point p2(1.2f, 0.0f, 1.7f);
  REQUIRE(!region.in(p2));
  Point p3(-2.1f, 0.05f, -1.0f);
  REQUIRE(!region.in(p3));

  std::mt19937 rng;

  for (int i = 0; i < 100; i++) {
    auto event = region.random_event(rng);
    auto direction = event.direction;

    REQUIRE(((direction.z <= std::sqrt(3) / 2) &&
             (direction.z >= -std::sqrt(3) / 2)));
  }
}

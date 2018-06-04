#include <random>

#include "util/test.h"
#include "event_generator.h"
#include "common/types.h"

TEST("3d/geometry/distribution/spherical_distribution") {
  using Distribution = PET3D::Distribution::SphericalDistribution<F>;
  using Vector = Distribution::Vector;

  std::mt19937 rng;

  Distribution distribution;
  for (int i = 0; i < 256; i++) {
    Vector dir = distribution(rng);

    CHECK(((std::abs(dir.x) <= 1) && (std::abs(dir.y) <= 1) &&
           (std::abs(dir.z) <= 1)));

    CHECK((dir.x * dir.x + dir.y * dir.y + dir.z * dir.z) == 1.0_e7);
  }
}

TEST("3d/geometry/distribution/cylinder_event_generator") {
  using Distribution = PET3D::Distribution::CylinderPointDistribution<F>;
  using Point = Distribution::Point;
  using F = Distribution::F;

  std::mt19937 rng;

  F radius = 2.0;
  F height = 3.0;
  Distribution distribution(radius, height);
  for (int i = 0; i < 100; i++) {
    Point p = distribution(rng);
    REQUIRE((p.x * p.x + p.y * p.y) <= radius * radius);
    REQUIRE(p.z <= height / 2);
    REQUIRE(p.z >= -height / 2);
  }
}

TEST("3d/geometry/distribution/ball_event_generator") {
  using Distribution = PET3D::Distribution::BallPointDistribution<F>;
  using Point = Distribution::Point;
  using F = Distribution::F;

  std::mt19937 rng;

  F radius = 2.0;
  Distribution distribution(radius);
  for (int i = 0; i < 100; i++) {
    Point p = distribution(rng);
    REQUIRE((p.x * p.x + p.y * p.y + p.z * p.z) <= radius * radius);
  }
}

TEST("3d/geometry/distribution/ellipsoid_event_generator") {
  using Distribution = PET3D::Distribution::EllipsoidPointDistribution<F>;
  using Point = Distribution::Point;
  using F = Distribution::F;

  std::mt19937 rng;

  F a = 3.0, b = 4.0, c = 5.0;
  Distribution distribution(a, b, c);
  for (int i = 0; i < 100; i++) {
    Point p = distribution(rng);
    REQUIRE((p.x * p.x / (a * a) + p.y * p.y / (b * b) + p.z * p.z / (c * c)) <=
            1);
  }
}

#include "util/test.h"

#include "event_generator.h"

using namespace PET3D;

TEST("3d/geometry/event_generator/spherical_distribution") {
  using Vector = SphericalDistribution<float>::Vector;

  SphericalDistribution<float> random_direction;
  std::mt19937 rng;

  for (int i = 0; i < 256; i++) {
    Vector dir = random_direction(rng);

    CHECK(((std::abs(dir.x) <= 1) && (std::abs(dir.y) <= 1) &&
           (std::abs(dir.z) <= 1)));

    CHECK((dir.x * dir.x + dir.y * dir.y + dir.z * dir.z) == 1.0_e7);
  }
}

TEST("3d/geometry/event_generator/voxel_event_generator") {

  using Event = VoxelEventGenerator<float>::Event;
  using Vector = VoxelEventGenerator<float>::Vector;
  using Point = VoxelEventGenerator<float>::Point;

  std::mt19937 rng;
  VoxelEventGenerator<float> event_generator(Point(1.0f, 2.0f, 3.0f),
                                             Vector(0.1f, 0.2f, 0.3f));

  for (int i = 0; i < 256; i++) {
    Event event = event_generator(rng);
    Point origin = event.origin;

    CHECK(((1.0 <= origin.x) && (origin.x <= 1.1)));
    CHECK(((2.0 <= origin.y) && (origin.y <= 2.2)));
    CHECK(((3.0 <= origin.z) && (origin.z <= 3.3)));

    Vector dir = event.direction;

    CHECK(((std::abs(dir.x) <= 1) && (std::abs(dir.y) <= 1) &&
           (std::abs(dir.z) <= 1)));

    CHECK((dir.x * dir.x + dir.y * dir.y + dir.z * dir.z) == 1.0_e7);
  }
}

TEST("3d/geometry/event_generator/cylinder_event_generator") {
  using Generator = PET3D::CylinderPointDistribution<float>;
  using Point = Generator::Point;
  using F = Generator::F;

  std::mt19937 rng;

  F radius = 2.0;
  F height = 3.0;
  Generator generator(radius, height);
  for (int i = 0; i < 100; i++) {
    Point p = generator(rng);
    REQUIRE((p.x * p.x + p.y * p.y) <= radius * radius);
    REQUIRE(p.z <= height / 2);
    REQUIRE(p.z >= -height / 2);
  }
}

TEST("3d/geometry/event_generator/ball_event_generator") {
  using Generator = PET3D::BallPointDistribution<float>;
  using Point = Generator::Point;
  using F = Generator::F;

  std::mt19937 rng;

  F radius = 2.0;
  Generator generator(radius);
  for (int i = 0; i < 100; i++) {
    Point p = generator(rng);
    REQUIRE((p.x * p.x + p.y * p.y + p.z * p.z) <= radius * radius);
  }
}

TEST("3d/geometry/event_generator/ellipsoid_event_generator") {
  using Generator = PET3D::EllipsoidPointDistribution<float>;
  using Point = Generator::Point;
  using F = Generator::F;

  std::mt19937 rng;

  F a = 3.0, b = 4.0, c = 5.0;
  Generator generator(a, b, c);
  for (int i = 0; i < 100; i++) {
    Point p = generator(rng);
    REQUIRE((p.x * p.x / (a * a) + p.y * p.y / (b * b) + p.z * p.z / (c * c)) <=
            1);
  }
}

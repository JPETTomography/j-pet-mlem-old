#include "util/test.h"

#include "event_generator.h"

using namespace PET3D;

TEST("3d/geometry/event_generator/spherical_distribution") {
  using Vector = spherical_distribution<float>::Vector;

  spherical_distribution<float> random_direction;
  std::mt19937 rng;

  for (int i = 0; i < 256; i++) {
    Vector dir = random_direction(rng);

    CHECK(((std::abs(dir.x) <= 1) && (std::abs(dir.y) <= 1) &&
           (std::abs(dir.z) <= 1)));

    CHECK((dir.x * dir.x + dir.y * dir.y + dir.z * dir.z) == 1.0_e7);
  }
}

TEST("3d/geometry/event_generator/voxel_event_generator") {

  using Event = voxel_event_generator<float>::Event;
  using Vector = voxel_event_generator<float>::Vector;
  using Point = voxel_event_generator<float>::Point;

  std::mt19937 rng;
  voxel_event_generator<float> event_generator(Point(1.0f, 2.0f, 3.0f),
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

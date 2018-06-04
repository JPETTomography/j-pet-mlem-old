#include <random>

#include "util/test.h"
#include "event_generator.h"
#include "common/types.h"

using VoxelEventGenerator = PET3D::VoxelEventGenerator<F>;
using Event = VoxelEventGenerator::Event;
using Vector = VoxelEventGenerator::Vector;
using Point = VoxelEventGenerator::Point;

TEST("3d/geometry/event_generator/voxel_event_generator") {
  std::mt19937 rng;
  VoxelEventGenerator event_generator(Point(1, 2, 3), Vector(0.1, 0.2, 0.3));

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

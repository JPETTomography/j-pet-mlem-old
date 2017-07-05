#include "util/test.h"

#include "2d/gate/gate_volume.h"
#include "2d/gate/gate_scanner_builder.h"
#include "2d/geometry/vector.h"
#include "common/types.h"

TEST("2d Gate volume") {
  using Box = Gate::D2::Box<F>;
  using Vector = Box::Vector;

  SECTION("Constructing") {
    Gate::D2::Box<float> world(1, 1);
    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(&world, &scanner);

    CHECK(scanner.size() == 0);
  }

  SECTION("One detector") {
    auto world = new Box(1, 1);
    auto detector = new Box(0.006, 0.024);
    detector->attach_crystal_sd();
    world->attach_daughter(detector);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    CHECK(scanner.size() == 1);
    auto d = scanner[0];
    CHECK(d.width() == Approx(0.006));
  }

  SECTION("One translated detector") {
    auto world = new Box(1, 1);
    auto detector = new Box(0.024, 0.006);
    detector->attach_crystal_sd();
    detector->set_translation(Vector(0.34, 0));
    world->attach_daughter(detector);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    CHECK(scanner.size() == 1);
    auto d = scanner[0];
    CHECK(d.width() == Approx(0.024));
    auto c = d.center();
    CHECK(c.x == Approx(0.34));
    CHECK(c.y == Approx(0.0));
  }

  SECTION("One translated && rotated detector") {
    auto world = new Box(1, 1);
    auto detector = new Box(0.006, 0.024);
    detector->attach_crystal_sd();
    detector->set_rotation(M_PI / 4);
    detector->set_translation(Vector(0.34, 0));
    world->attach_daughter(detector);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    CHECK(scanner.size() == 1);
    auto d = scanner[0];
    auto c = d.center();
    CHECK(c.x == Approx(0.34));
    CHECK(c.y == Approx(0.0));
  }
}

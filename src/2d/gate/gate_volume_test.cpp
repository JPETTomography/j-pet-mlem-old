#include "util/test.h"

#include "2d/gate/gate_volume.h"
#include "2d/gate/gate_scanner_builder.h"
#include "common/types.h"

TEST("2d Gate volume") {
  using Box = Gate::D2::Box<F>;

  SECTION("Constructing") {
    Gate::D2::Box<float> world;
    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(&world, &scanner);

    CHECK(scanner.size() == 0);
  }

  SECTION("One detector") {
    auto world = new Box();
    auto detector = new Box();
    detector->attach_crystal_sd();
    world->attach_daughter(detector);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    CHECK(scanner.size() == 1);
  }
}

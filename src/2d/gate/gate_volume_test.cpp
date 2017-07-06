#include <iostream>
#include <sstream>

#include "util/test.h"

#include "2d/gate/gate_volume.h"
#include "2d/gate/gate_scanner_builder.h"
#include "2d/geometry/vector.h"
#include "common/types.h"

#include "util/svg_ostream.h"

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

  SECTION("nested example") {
    auto world = new Box(1, 1);

    auto box1 = new Box(0.25, 0.25);
    box1->set_translation(Vector(0.15, 0));
    world->attach_daughter(box1);

    auto box2 = new Box(0.2, 0.1);
    box2->set_rotation(M_PI / 4);
    box1->attach_daughter(box2);

    auto da = new Box(0.05, 0.1);
    da->set_translation(Vector(0.05, 0));
    da->attach_crystal_sd();
    auto db = new Box(0.05, 0.1);
    db->set_translation(Vector(-0.05, 0));
    db->attach_crystal_sd();

    box2->attach_daughter(da);
    box2->attach_daughter(db);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    CHECK(scanner.size() == 2);

    util::svg_ostream<F> out(
        "gate_volume_scanner_test.svg", 1., 1., 1024., 1024l);
    out << scanner;

    std::ifstream test_in("src/2d/geometry/gate_volume_test.txt");
    if (test_in) {
      auto d_a = scanner[0];
      for (int i = 0; i < 4; i++) {
        F x, y;
        test_in >> x >> y;
        CHECK(d_a[i].x == Approx(x));
        CHECK(d_a[i].y == Approx(y));
      }

      auto d_b = scanner[1];
      for (int i = 0; i < 4; i++) {
        F x, y;
        test_in >> x >> y;
        CHECK(d_b[i].x == Approx(x));
        CHECK(d_b[i].y == Approx(y));
      }

    } else {
      WARN("could not open file `src/2d/geometry/gate_volume_test.txt'");
    }
  }

  SECTION("Simple linear repeater") {
    auto world = new Box(1, 1);
    auto da = new Box(0.05, 0.1);

    da->attach_repeater(new Gate::D2::Linear<F>(4, Vector(0.1, 0.0)));
    da->attach_crystal_sd();

    world->attach_daughter(da);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    util::svg_ostream<F> out(
        "gate_volume_s_repeater_test.svg", 1., 1., 1024., 1024l);
    out << scanner;
    REQUIRE(scanner.size() == 4);

    std::ifstream test_in("src/2d/geometry/gate_volume_s_repeater_test.txt");
    if (test_in) {
      for (int j = 0; j < 4; j++) {
        auto d = scanner[j];
        for (int i = 0; i < 4; i++) {
          F x, y;
          test_in >> x >> y;
          REQUIRE(d[i].x == Approx(x));
          REQUIRE(d[i].y == Approx(y));
        }
      }

    } else {
      WARN(
          "could not open file "
          "`src/2d/geometry/gate_volume_s_repeater_test.txt'");
    }
  }
}

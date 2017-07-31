#include <iostream>
#include <sstream>

#include "util/test.h"

#include "2d/gate/gate_volume.h"
#include "2d/gate/gate_scanner_builder.h"
#include "2d/geometry/vector.h"
#include "common/types.h"
#include "common/mathematica_graphics.h"
#include "2d/barrel/scanner_builder.h"

#include "gate_volume_builder.h"

#include "util/svg_ostream.h"

TEST("2d Gate volume") {
  using Box = Gate::D2::Box<F>;
  using Vector = Box::Vector;
  using Cylinder = Gate::D2::Cylinder<F>;

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

  SECTION("Simple ring repeater") {
    auto world = new Box(1, 1);
    auto da = new Box(0.05, 0.1);

    da->attach_repeater(new Gate::D2::Ring<F>(5, Vector(0.0, 0.0)));
    da->attach_crystal_sd();
    da->set_translation(Vector(0, 0.2));

    world->attach_daughter(da);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    util::svg_ostream<F> out(
        "gate_volume_s_repeater_ring_test.svg", 1., 1., 1024., 1024l);
    out << scanner;
    REQUIRE(scanner.size() == 5);

    std::ifstream test_in(
        "src/2d/geometry/gate_volume_s_repeater_ring_test.txt");
    if (test_in) {
      for (int j = 0; j < 5; j++) {
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
          "`src/2d/geometry/gate_volume_s_repeater_ring_test.txt'");
    }
  }

  SECTION("Simple ring repeater off center") {
    auto world = new Box(1, 1);
    auto mod = new Box(1, 1);

    auto da = new Box(0.05, 0.1);
    world->attach_daughter(mod);

    da->attach_repeater(new Gate::D2::Ring<F>(5, Vector(0.2, 0.10)));
    da->attach_crystal_sd();
    da->set_translation(Vector(0, 0.2));

    mod->attach_daughter(da);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner;
    builder.build(world, &scanner);

    util::svg_ostream<F> out(
        "gate_volume_s_repeater_ring_off_test.svg", 1., 1., 1024., 1024l);
    out << scanner;
    REQUIRE(scanner.size() == 5);

    std::ifstream test_in(
        "src/2d/geometry/gate_volume_s_repeater_ring_off_test.txt");
    if (test_in) {
      for (int j = 0; j < 5; j++) {
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
          "`src/2d/geometry/gate_volume_s_repeater_ring_off_test.txt'");
    }
  }

  SECTION("new modules") {
    auto world = new Box(2, 2);

    auto layer_new = new Cylinder(0.35, 0.4);
    world->attach_daughter(layer_new);

    auto module = new Box(0.026, 0.0085);
    module->set_translation(Vector(0.37236, 0));
    module->attach_repeater(new Gate::D2::Ring<F>(24, Vector(0.0, 0.0)));

    layer_new->attach_daughter(module);

    auto scintillator = new Box(0.024, 0.006);

    scintillator->attach_repeater(
        new Gate::D2::Linear<F>(13, Vector(0, 0.007)));
    scintillator->attach_crystal_sd();

    module->attach_daughter(scintillator);

    Gate::D2::GenericScannerBuilder<F, S> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S> scanner(
        0.4, 0.8);
    builder.build(world, &scanner);
    REQUIRE(13 * 24 == scanner.size());

    util::svg_ostream<F> out(
        "gate_volume_new_modules.svg", .9, .9, 1024., 1024l);
    out << scanner;

    std::ofstream mout("gate_volume_new_modules.m");
    Common::MathematicaGraphics<F> mgraphics(mout);
    mgraphics.add(scanner);
  }

  SECTION("new full") {
    auto world = new Box(2, 2);

    auto layer_new = new Cylinder(0.35, 0.4);
    world->attach_daughter(layer_new);

    auto module = new Box(0.026, 0.0085);
    module->set_translation(Vector(0.37236, 0));
    module->attach_repeater(new Gate::D2::Ring<F>(24, Vector(0.0, 0.0)));

    layer_new->attach_daughter(module);

    auto scintillator = new Box(0.024, 0.006);
    scintillator->attach_repeater(
        new Gate::D2::Linear<F>(13, Vector(0, 0.007)));
    scintillator->attach_crystal_sd();

    module->attach_daughter(scintillator);

    auto layer_1 = new Cylinder(0.425 - 0.005, 0.425 + 0.005);
    world->attach_daughter(layer_1);
    auto scintillator_1 = new Box(0.021, 0.009);
    scintillator_1->set_translation(Vector(0.425, 0));
    scintillator_1->attach_repeater(new Gate::D2::Ring<F>(48, Vector(0, 0)));
    scintillator_1->attach_crystal_sd();
    layer_1->attach_daughter(scintillator_1);

    auto layer_2 = new Cylinder(0.4675 - 0.005, 0.4675 + 0.005);
    world->attach_daughter(layer_2);
    auto scintillator_2 = new Box(0.021, 0.009);
    scintillator_2->set_translation(Vector(0.4675, 0));
    scintillator_2->attach_repeater(
        new Gate::D2::Ring<F>(48, Vector(0, 0), M_PI / 48));
    scintillator_2->attach_crystal_sd();
    layer_2->attach_daughter(scintillator_2);

    auto layer_3 = new Cylinder(0.575 - 0.005, 0.575 + 0.005);
    world->attach_daughter(layer_3);
    auto scintillator_3 = new Box(0.021, 0.009);
    scintillator_3->set_translation(Vector(0.575, 0));
    scintillator_3->attach_repeater(
        new Gate::D2::Ring<F>(96, Vector(0, 0), M_PI / 96));
    scintillator_3->attach_crystal_sd();
    layer_3->attach_daughter(scintillator_3);

    Gate::D2::GenericScannerBuilder<F, S, 512> builder;
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S, 512>
        scanner(0.4, 0.8);
    builder.build(world, &scanner);
    REQUIRE(13 * 24 + 48 + 48 + 96 == scanner.size());

    util::svg_ostream<F> out("new_full.svg", .9, .9, 1024., 1024l);
    out << scanner;

    std::ofstream mout("new_full.m");
    Common::MathematicaGraphics<F> mgraphics(mout);
    mgraphics.add(scanner);

    auto n_detectors = builder.count_cristals(world);
    REQUIRE(n_detectors == scanner.size());
  }

  SECTION("new full from builder") {

    auto world = Gate::D2::build_new_full_scanner_volume<F>();

    Gate::D2::GenericScannerBuilder<F, S, 512> builder;
    auto scanner = builder.build_with_8_symmetries(world);
    REQUIRE(13 * 24 + 48 + 48 + 96 == scanner.size());

    util::svg_ostream<F> out("new_full_from_builder.svg", .9, .9, 1024., 1024l);
    out << scanner;

    std::ofstream mout("new_full_from_builder.m");
    Common::MathematicaGraphics<F> mgraphics(mout);
    mgraphics.add(scanner);

    auto n_detectors = builder.count_cristals(world);
    REQUIRE(n_detectors == scanner.size());
  }
}

TEST("old multi ring") {
  using Box = Gate::D2::Box<F>;
  using Vector = Box::Vector;
  using Cylinder = Gate::D2::Cylinder<F>;

  Gate::D2::Volume<F>* world = new Box(2, 2);

  auto layer_1 = new Cylinder(0.425 - 0.005, 0.425 + 0.005);
  world->attach_daughter(layer_1);
  auto scintillator_1 = new Box(0.021, 0.009);
  scintillator_1->set_translation(Vector(0.425, 0));
  scintillator_1->attach_repeater(new Gate::D2::Ring<F>(48, Vector(0, 0)));
  scintillator_1->attach_crystal_sd();
  layer_1->attach_daughter(scintillator_1);

  auto layer_2 = new Cylinder(0.4675 - 0.005, 0.4675 + 0.005);
  world->attach_daughter(layer_2);
  auto scintillator_2 = new Box(0.021, 0.009);
  scintillator_2->set_translation(Vector(0.4675, 0));
  scintillator_2->attach_repeater(
      new Gate::D2::Ring<F>(48, Vector(0, 0), M_PI / 48));
  scintillator_2->attach_crystal_sd();
  layer_2->attach_daughter(scintillator_2);

  auto layer_3 = new Cylinder(0.575 - 0.005, 0.575 + 0.005);
  world->attach_daughter(layer_3);
  auto scintillator_3 = new Box(0.021, 0.009);
  scintillator_3->set_translation(Vector(0.575, 0));
  scintillator_3->attach_repeater(
      new Gate::D2::Ring<F>(96, Vector(0, 0), M_PI / 96));
  scintillator_3->attach_crystal_sd();
  layer_3->attach_daughter(scintillator_3);

  Gate::D2::GenericScannerBuilder<F, S, 512> builder;

  world = Gate::D2::build_big_barrel_volume<F>();

  SECTION("build") {

    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S, 512>
        scanner(0.4, 0.8);
    builder.build(world, &scanner);
    REQUIRE(48 + 48 + 96 == scanner.size());

    util::svg_ostream<F> out("new_full.svg", .9, .9, 1024., 1024l);
    out << scanner;

    std::ofstream mout("new_full.m");
    Common::MathematicaGraphics<F> mgraphics(mout);
    mgraphics.add(scanner);
  }

  SECTION("build_with*_symmetries") {
    auto scanner = builder.build_with_8_symmetries(world);
    REQUIRE(48 + 48 + 96 == scanner.size());

    auto s_descriptor = scanner.symmetry_descriptor();

    std::vector<F> rads = { 1, 2, 3 };
    std::vector<F> rots = { 0, 0.5, 0.5 };
    std::vector<int> n_dets = { 48, 48, 96 };

    auto ref_scanner = PET2D::Barrel::ScannerBuilder<
        PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>,
                                      S,
                                      512>>::build_multiple_rings(rads,
                                                                  rots,
                                                                  n_dets,
                                                                  0.009,
                                                                  0.021);

    for (S s = 0; s < 8; s++) {
      for (S d = 0; d < scanner.size(); d++) {
        REQUIRE(s_descriptor.symmetric_detector(d, s) ==
                ref_scanner.symmetry_descriptor().symmetric_detector(d, s));
      }
    }
  }
}

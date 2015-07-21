#include <cmath>

#include "util/test.h"

#include "3d/hybrid/scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/square_detector.h"
#include "3d/geometry/point.h"
#include "3d/geometry/vector.h"

#include "common/model.h"
#include "common/types.h"

F center_radius = 0.180;
F scintillator_height = 0.019;
F scintillator_width = 0.05;
F inner_radius = center_radius - scintillator_height / 2;
F length = 0.30;

F minimal_angle = std::atan2(inner_radius, length / 2);

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, S, 24>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

TEST("3d/hybrid/detector_set/escape_through_endcap") {
  Scanner2D scanner_2d(inner_radius, scintillator_height);
  Scanner scanner(scanner_2d, length);

  {
    PET3D::Event<F> event(Point(0, 0, 0), Vector(0, 0, 1));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, M_PI / 2));
    CHECK(!scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(1, M_PI / 2));
    CHECK(!scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, M_PI / 4));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, 0.99f * minimal_angle));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, 1.01f * minimal_angle));
    CHECK(!scanner.escapes_through_endcap(event));
  }
}

TEST("3d/hybrid/detector_set/detect", "detect") {
  Scanner2D scanner_2d =
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
          inner_radius, 24, scintillator_height, scintillator_width);
  Scanner scanner(scanner_2d, length);
  Common::AlwaysAccept<F> model;

  {
    PET3D::Event<F> event(Point(0, 0, 0), Vector(0, 0, 1));

    Scanner::Response response;

    CHECK(!scanner.detect(model, model, event, response));
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, M_PI / 2));

    Scanner::Response response;

    REQUIRE(scanner.detect(model, model, event, response));
    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
    CHECK(response.z_up == 0.0_e7);
    CHECK(response.z_dn == 0.0_e7);
    CHECK(response.dl == 0.0_e7);
  }

  {
    PET3D::Event<F> event(Point(0, 0, 0),
                          Vector::from_euler_angles(0, M_PI / 2.5f));

    Scanner::Response response;

    CHECK(scanner.detect(model, model, event, response));

    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
    F z = inner_radius / std::tan(M_PI / 2.5);
    CHECK(response.z_up == Approx(-z).epsilon(1e-7));
    CHECK(response.z_dn == Approx(z).epsilon(1e-7));
    CHECK(response.dl == 0.0_e7);
  }
}

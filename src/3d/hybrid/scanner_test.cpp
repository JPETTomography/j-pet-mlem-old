#include <cmath>

#include "util/test.h"

#include "3d/hybrid/scanner.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/square_detector.h"
#include "3d/geometry/point.h"
#include "3d/geometry/vector.h"
#include "common/model.h"

float center_radius = 0.180f;
float scintillator_height = 0.019f;
float scintillator_width = 0.05f;
float inner_radius = center_radius - scintillator_height / 2;
float length = 0.30;

float minimal_angle = std::atan2(inner_radius, length / 2);

using SquareDetector = PET2D::Barrel::SquareDetector<float>;
using Scanner2D = PET2D::Barrel::GenericScanner<SquareDetector, short, 24>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Vector = PET3D::Vector<float>;
using Point = PET3D::Point<float>;

TEST("3d/hybrid/detector_set/escape_through_endcap") {
  Scanner2D scanner_2d(inner_radius, scintillator_height);
  Scanner scanner(scanner_2d, length);

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector(0.0f, 0.0f, 1.0f));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.0));
    CHECK(!scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(1.0f, M_PI / 2.0));
    CHECK(!scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 4.0));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(
        Point(0.0f, 0.0f, 0.0f),
        Vector::from_euler_angles(0.0f, 0.99f * minimal_angle));
    CHECK(scanner.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(
        Point(0.0f, 0.0f, 0.0f),
        Vector::from_euler_angles(0.0f, 1.01f * minimal_angle));
    CHECK(!scanner.escapes_through_endcap(event));
  }
}

TEST("3d/hybrid/detector_set/detect", "detect") {
  Scanner2D scanner_2d =
      PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
          inner_radius, 24, scintillator_height, scintillator_width);
  Scanner scanner(scanner_2d, length);
  Common::AlwaysAccept<float> model;

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector(0.0f, 0.0f, 1.0f));

    Scanner::Response response;

    CHECK(!scanner.detect(model, model, event, response));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.0));

    Scanner::Response response;

    REQUIRE(scanner.detect(model, model, event, response));
    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
#if DEBUG
    std::cerr << lor.first << " " << lor.second << "\n";
#endif
    CHECK(response.z_up == 0.0_e7);
    CHECK(response.z_dn == 0.0_e7);
    CHECK(response.dl == 0.0_e7);
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.5f));

    Scanner::Response response;

    CHECK(scanner.detect(model, model, event, response));

    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
#if DEBUG
    std::cerr << lor.first << " " << lor.second << "\n";
#endif
    float z = inner_radius / std::tan(M_PI / 2.5);
#if DEBUG
    std::cerr << z << "\n";
#endif
    CHECK(response.z_up == Approx(-z).epsilon(1e-7));
    CHECK(response.z_dn == Approx(z).epsilon(1e-7));
    CHECK(response.dl == 0.0_e7);
  }
}

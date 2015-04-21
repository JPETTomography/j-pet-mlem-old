#include <cmath>

#include "util/test.h"

#include "3d/longitudinal/detector_set.h"
#include "2d/barrel/detector_set.h"
#include "2d/barrel/detector_set_builder.h"
#include "2d/barrel/square_detector.h"
#include "3d/geometry/point.h"
#include "3d/geometry/vector.h"
#include "2d/barrel/model.h"

float center_radius = 0.180f;
float scintillator_height = 0.019f;
float scintillator_width = 0.05f;
float inner_radius = center_radius - scintillator_height / 2;
float length = 0.30;

float minimal_angle = std::atan2(inner_radius, length / 2);

TEST("3d/longitudinal/detector_set/escape_through_endcap") {
  using SquareDetector = PET2D::Barrel::SquareDetector<float>;
  using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareDetector,24>;
  using DetectorSet = PET3D::Longitudinal::DetectorSet<DetectorSet2D>;
  using Vector = PET3D::Vector<float>;
  using Point = PET3D::Point<float>;

  DetectorSet2D detector_set_2D(inner_radius, scintillator_height);
  DetectorSet detector_set(detector_set_2D, length);

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector(0.0f, 0.0f, 1.0f));
    CHECK(detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.0));
    CHECK(!detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(1.0f, M_PI / 2.0));
    CHECK(!detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 4.0));
    CHECK(detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(
        Point(0.0f, 0.0f, 0.0f),
        Vector::from_euler_angles(0.0f, 0.99f * minimal_angle));
    CHECK(detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(
        Point(0.0f, 0.0f, 0.0f),
        Vector::from_euler_angles(0.0f, 1.01f * minimal_angle));
    CHECK(!detector_set.escapes_through_endcap(event));
  }
}

TEST("3d/longitudinal/detector_set/detect", "detect") {
  using SquareDetector = PET2D::Barrel::SquareDetector<float>;
  using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareDetector, 24>;
  using DetectorSet = PET3D::Longitudinal::DetectorSet<DetectorSet2D>;
  using Vector = PET3D::Vector<float>;
  using Point = PET3D::Point<float>;

  DetectorSet2D detector_set_2D =
      PET2D::Barrel::DetectorSetBuilder<DetectorSet2D>::buildSingleRing(
          inner_radius, 24, scintillator_height, scintillator_width);
  DetectorSet detector_set(detector_set_2D, length);
  PET2D::Barrel::AlwaysAccept<> model;

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector(0.0f, 0.0f, 1.0f));

    DetectorSet::Response response;

    CHECK(!detector_set.detect(model, model, event, response));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.0));

    DetectorSet::Response response;

    REQUIRE(detector_set.detect(model, model, event, response));
    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
    // std::cerr << lor.first << " " << lor.second << "\n";
    CHECK(response.z_up == 0.0_e7);
    CHECK(response.z_dn == 0.0_e7);
    CHECK(response.dl == 0.0_e7);
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, M_PI / 2.5f));

    DetectorSet::Response response;

    CHECK(detector_set.detect(model, model, event, response));

    auto lor = response.lor;
    CHECK(lor.first == 12);
    CHECK(lor.second == 0);
    // std::cerr << lor.first << " " << lor.second << "\n";
    float z = inner_radius / std::tan(M_PI / 2.5);
    //    std::cerr<<z<<"\n";
    CHECK(response.z_up == Approx(-z).epsilon(1e-7));
    CHECK(response.z_dn == Approx(z).epsilon(1e-7));
    CHECK(response.dl == 0.0_e7);
  }
}

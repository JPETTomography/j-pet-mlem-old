#include <cmath>

#include "util/test.h"

#include "3d/longitudinal/detector_set.h"
#include "2d/barrel/detector_set.h"
#include "2d/barrel/square_detector.h"
#include "3d/geometry/point.h"
#include "3d/geometry/vector.h"
#include "2d/barrel/model.h"


TEST("3d/longitudinal/detector_set/escape_through_endcap") {
  using SquareDetector = PET2D::Barrel::SquareDetector<float>;
  using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareDetector>;
  using DetectorSet = PET3D::Longitudinal::DetectorSet<DetectorSet2D, 24>;
  using Vector = PET3D::Vector<float>;
  using Point = PET3D::Point<float>;

  DetectorSet2D detector_set_2D(0.17f, 0.019f);

  DetectorSet detector_set(detector_set_2D, 0.5f);

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
    CHECK(!detector_set.escapes_through_endcap(event));
  }

  {
    PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f),
                              Vector::from_euler_angles(0.0f, 0.59f));
    CHECK(detector_set.escapes_through_endcap(event));
  }
}

TEST("3d/longitudinal/detector_set/detect", "detect") {
  using SquareDetector = PET2D::Barrel::SquareDetector<float>;
  using DetectorSet2D = PET2D::Barrel::DetectorSet<SquareDetector>;
  using DetectorSet = PET3D::Longitudinal::DetectorSet<DetectorSet2D, 24>;
  using Vector = PET3D::Vector<float>;
  using Point = PET3D::Point<float>;

  DetectorSet2D detector_set_2D(0.17f, 0.019f);

  DetectorSet detector_set(detector_set_2D, 0.5f);
  PET2D::Barrel::AlwaysAccept<> model;
  PET3D::Event<float> event(Point(0.0f, 0.0f, 0.0f), Vector(0.0f, 0.0f, 1.0f));
  PET2D::Barrel::LOR<> lor;
  float position;
  if(!detector_set.escapes_through_endcap(event)) {

  }

}

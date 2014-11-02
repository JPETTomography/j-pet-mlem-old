#include <cmath>

#include "util/test.h"

#include "model.h"
#include "compound_detector.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/compound_detector/math") {
  SECTION("square_detector") {
    CompoundDetector<SquareDetector<>> detector;
    detector.emplace_back(1., 1., 2., 2.);
    detector.emplace_back(1., 5., 2., 2.);
    detector.emplace_back(5., 1., 2., 2.);
    detector.emplace_back(5., 5., 2., 2.);

    CHECK(Point<>(1., 1.) == detector.circumscribed(0).center);
    CHECK(Point<>(1., 5.) == detector.circumscribed(1).center);
    CHECK(Point<>(5., 1.) == detector.circumscribed(2).center);
    CHECK(Point<>(5., 5.) == detector.circumscribed(3).center);

    CHECK(std::sqrt(2.) == detector.circumscribed(0).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(1).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(2).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(3).radius);
  }

  SECTION("circle_detector") {
    CompoundDetector<CircleDetector<>> detector;
    detector.emplace_back(2., Point<>(1., 1.));
    detector.emplace_back(2., Point<>(1., 5.));
    detector.emplace_back(2., Point<>(5., 1.));
    detector.emplace_back(2., Point<>(5., 5.));

    CHECK(Point<>(1., 1.) == detector.circumscribed(0).center);
    CHECK(Point<>(1., 5.) == detector.circumscribed(1).center);
    CHECK(Point<>(5., 1.) == detector.circumscribed(2).center);
    CHECK(Point<>(5., 5.) == detector.circumscribed(3).center);

    CHECK(2. == detector.circumscribed(0).radius);
    CHECK(2. == detector.circumscribed(1).radius);
    CHECK(2. == detector.circumscribed(2).radius);
    CHECK(2. == detector.circumscribed(3).radius);
  }
}

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

    CHECK(Point<>(1., 1.) == detector.circumscribed(0));
    CHECK(Point<>(1., 5.) == detector.circumscribed(1));
    CHECK(Point<>(5., 1.) == detector.circumscribed(2));
    CHECK(Point<>(5., 5.) == detector.circumscribed(3));

    CHECK(std::sqrt(2.) == detector.circumscribed(0).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(1).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(2).radius);
    CHECK(std::sqrt(2.) == detector.circumscribed(3).radius);

    LOR<> lor;
    double position;
    Event<> e1(0, 0, M_PI_2);
    CHECK(0 == detector.detect(e1, e1, e1, lor, position));

    Event<> e2(4, 4, M_PI_2);
    CHECK(2 == detector.detect(e2, e2, e2, lor, position));
  }

  SECTION("circle_detector") {
    CompoundDetector<CircleDetector<>> detector;
    detector.emplace_back(2., Point<>(1., 1.));
    detector.emplace_back(2., Point<>(1., 5.));
    detector.emplace_back(2., Point<>(5., 1.));
    detector.emplace_back(2., Point<>(5., 5.));

    CHECK(Point<>(1., 1.) == detector.circumscribed(0));
    CHECK(Point<>(1., 5.) == detector.circumscribed(1));
    CHECK(Point<>(5., 1.) == detector.circumscribed(2));
    CHECK(Point<>(5., 5.) == detector.circumscribed(3));

    CHECK(2. == detector.circumscribed(0).radius);
    CHECK(2. == detector.circumscribed(1).radius);
    CHECK(2. == detector.circumscribed(2).radius);
    CHECK(2. == detector.circumscribed(3).radius);
  }
}

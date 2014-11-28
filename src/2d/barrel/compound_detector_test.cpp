#include <cmath>

#include "util/test.h"

#include "model.h"
#include "compound_detector.h"
#include "detector_ring.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/compound_detector/math") {
  SECTION("square_detector") {
    CompoundDetector<SquareDetector<>> detector;
    detector.emplace_back(1., 1., 2., 2.);  // 0
    detector.emplace_back(1., 5., 2., 2.);  // 1
    detector.emplace_back(5., 1., 2., 2.);  // 2
    detector.emplace_back(5., 5., 2., 2.);  // 3

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

    SECTION("horizontal") {
      Event<> e1(3, 1, 0);
      CHECK(0 == detector.close_indices(e1).first[0]);
      CHECK(2 == detector.close_indices(e1).second[0]);
      Event<> e2(3, 4, 0);
      CHECK(1 == detector.close_indices(e2).first[0]);
      CHECK(3 == detector.close_indices(e2).second[0]);
    }

    SECTION("vertical") {
      Event<> e1(1, 3, M_PI_2);
      CHECK(0 == detector.close_indices(e1).first[0]);
      CHECK(1 == detector.close_indices(e1).second[0]);
      Event<> e2(4, 3, M_PI_2);
      CHECK(2 == detector.close_indices(e2).first[0]);
      CHECK(3 == detector.close_indices(e2).second[0]);
    }
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

TEST("2d/barrel/compound_detector/detect") {
  SECTION("two_rings") {
    DetectorRing<SquareDetector<>> inner_ring(16, 1., .1, .1);
    DetectorRing<SquareDetector<>> outer_ring(16, 1.4, .1, .1);
    CompoundDetector<SquareDetector<>> detector;
    for (auto& square_detector : inner_ring) {
      detector.push_back(square_detector);
    }
    for (auto& square_detector : outer_ring) {
      detector.push_back(square_detector);
    }
    CHECK(32 == detector.size());

    DetectorRing<>::LOR lor;
    AlwaysAccept<> model;
    double position;

    SECTION("horizontal") {
      Event<> e1(0, 0, 0);
      auto i1 = detector.close_indices(e1);
      CHECK(2 == i1.first.size());
      CHECK(2 == i1.second.size());
      CHECK(8 == i1.first[0]);
      CHECK(0 == i1.second[0]);
      CHECK(24 == i1.first[1]);
      CHECK(16 == i1.second[1]);

      CHECK(2 == detector.detect(model, model, e1, lor, position));

      Event<> e2(0, .050001, 0);
      auto i2 = detector.close_indices(e1);
      CHECK(2 == i2.first.size());
      CHECK(2 == i2.second.size());
      CHECK(8 == i2.first[0]);
      CHECK(0 == i2.second[0]);
      CHECK(24 == i2.first[1]);
      CHECK(16 == i2.second[1]);

      CHECK(0 == detector.detect(model, model, e2, lor, position));

      Event<> e3(0, .049999, 0);
      CHECK(2 == detector.detect(model, model, e3, lor, position));
    }
  }
}

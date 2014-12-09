#include <cmath>

#include "util/test.h"

#include "model.h"
#include "detector_set.h"
#include "detector_ring.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/detector_set/math") {
  SECTION("square_detector") {
    DetectorSet<SquareDetector<>> detector;
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

    SECTION("horizontal") {
      {
        Event<> e(3, 1, 0);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(0 == left[0]);
        CHECK(2 == right[0]);
      }
      {
        Event<> e(3, 4, 0);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(1 == left[0]);
        CHECK(3 == right[0]);
      }
    }

    SECTION("vertical") {
      {
        Event<> e(1, 3, M_PI_2);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(0 == left[0]);
        CHECK(1 == right[0]);
      }
      {
        Event<> e(4, 3, M_PI_2);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(2 == left[0]);
        CHECK(3 == right[0]);
      }
    }
  }

  SECTION("circle_detector") {
    DetectorSet<CircleDetector<>> detector;
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

TEST("2d/barrel/detector_set/detect") {
  SECTION("two_rings") {
    DetectorRing<SquareDetector<>> inner_ring(16, 1., .1, .1);
    DetectorRing<SquareDetector<>> outer_ring(16, 1.4, .1, .1);
    DetectorSet<SquareDetector<>> detector;
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
      {
        Event<> e(0, 0, 0);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(2 == left.size());
        CHECK(2 == right.size());
        CHECK(8 == left[0]);
        CHECK(0 == right[0]);
        CHECK(24 == left[1]);
        CHECK(16 == right[1]);

        CHECK(2 == detector.detect(model, model, e, lor, position));
      }
      {
        Event<> e(0, .050001, 0);
        DetectorSet<SquareDetector<>>::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(2 == left.size());
        CHECK(2 == right.size());
        CHECK(8 == left[0]);
        CHECK(0 == right[0]);
        CHECK(24 == left[1]);
        CHECK(16 == right[1]);

        CHECK(0 == detector.detect(model, model, e, lor, position));
      }
      {
        Event<> e(0, .049999, 0);
        CHECK(2 == detector.detect(model, model, e, lor, position));
      }
    }
  }
}

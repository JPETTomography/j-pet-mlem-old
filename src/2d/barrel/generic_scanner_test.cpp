#include <cmath>

#include "util/test.h"

#include "common/model.h"
#include "generic_scanner.h"
#include "ring_scanner.h"
#include "scanner_builder.h"

TEST("2d/barrel/detector_set/math") {
  SECTION("square_detector") {
    using SquareDetector = PET2D::Barrel::SquareDetector<double>;
    using Scanner = PET2D::Barrel::GenericScanner<SquareDetector, short, 24>;
    using Event = Scanner::Event;
    using Point = Scanner::Point;
    Scanner scanner;

    scanner.emplace_back(1., 1., 2., 2.);  // 0
    scanner.emplace_back(1., 5., 2., 2.);  // 1
    scanner.emplace_back(5., 1., 2., 2.);  // 2
    scanner.emplace_back(5., 5., 2., 2.);  // 3

    CHECK(Point(1., 1.) == scanner.circumscribed(0).center);
    CHECK(Point(1., 5.) == scanner.circumscribed(1).center);
    CHECK(Point(5., 1.) == scanner.circumscribed(2).center);
    CHECK(Point(5., 5.) == scanner.circumscribed(3).center);

    CHECK(std::sqrt(2.) == scanner.circumscribed(0).radius);
    CHECK(std::sqrt(2.) == scanner.circumscribed(1).radius);
    CHECK(std::sqrt(2.) == scanner.circumscribed(2).radius);
    CHECK(std::sqrt(2.) == scanner.circumscribed(3).radius);

    SECTION("horizontal") {
      {
        Event e(3, 1, 0);
        Scanner::Indices left, right;
        scanner.close_indices(e, left, right);
        CHECK(0 == left[0]);
        CHECK(2 == right[0]);
      }
      {
        Event e(3, 4, 0);
        Scanner::Indices left, right;
        scanner.close_indices(e, left, right);
        CHECK(1 == left[0]);
        CHECK(3 == right[0]);
      }
    }

    SECTION("vertical") {
      {
        Event e(1, 3, M_PI_2);
        Scanner::Indices left, right;
        scanner.close_indices(e, left, right);
        CHECK(0 == left[0]);
        CHECK(1 == right[0]);
      }
      {
        Event e(4, 3, M_PI_2);
        Scanner::Indices left, right;
        scanner.close_indices(e, left, right);
        CHECK(2 == left[0]);
        CHECK(3 == right[0]);
      }
    }
  }

  SECTION("circle_detector") {
    using CircleDetector = PET2D::Barrel::CircleDetector<float>;
    using Detector = PET2D::Barrel::GenericScanner<CircleDetector, short, 24>;
    using Point = Detector::Point;
    Detector detector;

    detector.emplace_back(2., Point(1., 1.));
    detector.emplace_back(2., Point(1., 5.));
    detector.emplace_back(2., Point(5., 1.));
    detector.emplace_back(2., Point(5., 5.));

    CHECK(Point(1., 1.) == detector.circumscribed(0).center);
    CHECK(Point(1., 5.) == detector.circumscribed(1).center);
    CHECK(Point(5., 1.) == detector.circumscribed(2).center);
    CHECK(Point(5., 5.) == detector.circumscribed(3).center);

    CHECK(2. == detector.circumscribed(0).radius);
    CHECK(2. == detector.circumscribed(1).radius);
    CHECK(2. == detector.circumscribed(2).radius);
    CHECK(2. == detector.circumscribed(3).radius);
  }
}

TEST("2d/barrel/detector_set/detect") {
  SECTION("two_rings") {

    using SquareDetector = PET2D::Barrel::SquareDetector<float>;
    using Detector = PET2D::Barrel::GenericScanner<SquareDetector, short, 128>;
    using Event = Detector::Event;

    using Response = typename Detector::Response;

    Detector inner_ring =
        PET2D::Barrel::ScannerBuilder<Detector>::build_single_ring(
            1., 16, .1, .1);
    Detector outer_ring =
        PET2D::Barrel::ScannerBuilder<Detector>::build_single_ring(
            1.4, 16, .1, .1);

    Detector detector;

    for (auto& square_detector : inner_ring) {
      detector.push_back(square_detector);
    }
    for (auto& square_detector : outer_ring) {
      detector.push_back(square_detector);
    }
    CHECK(32 == detector.size());

    Common::AlwaysAccept<float> model;
    Response response;
    SECTION("horizontal") {
      {
        Event e(0, 0, 0);
        Detector::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(2 == left.size());
        CHECK(2 == right.size());
        CHECK(8 == left[0]);
        CHECK(0 == right[0]);
        CHECK(24 == left[1]);
        CHECK(16 == right[1]);

        CHECK(2 == detector.detect(model, model, e, response));
      }
      {
        Event e(0, .050001, 0);
        Detector::Indices left, right;
        detector.close_indices(e, left, right);
        CHECK(2 == left.size());
        CHECK(2 == right.size());
        CHECK(8 == left[0]);
        CHECK(0 == right[0]);
        CHECK(24 == left[1]);
        CHECK(16 == right[1]);

        CHECK(0 == detector.detect(model, model, e, response));
      }
      {
        Event e(0, .049999, 0);
        CHECK(2 == detector.detect(model, model, e, response));
      }
    }
  }
}

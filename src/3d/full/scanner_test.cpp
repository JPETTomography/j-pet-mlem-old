#include <random>

#include "util/test.h"

#include "3d/full/scanner.h"
#include "common/model.h"

#include "common/types.h"

using Vector = PET3D::Vector<F>;
using Point = PET3D::Point<F>;

using RNG = std::mt19937;

TEST_CASE("3d/full/scanner") {
  PET3D::Full::Scanner<F, S> scanner;

  F length = 0.500;
  F width = 0.007;
  F height = 0.019;
  F R = 0.500;

  scanner.add_volume(PET3D::Full::BoxVolume<F>::AAB(
      Point(R - height / 2, -width / 2, -length / 2),
      Point(R + height / 2, width / 2, length / 2)));

  scanner.add_volume(PET3D::Full::BoxVolume<F>::AAB(
      Point(-R + height / 2, -width / 2, -length / 2),
      Point(-R - height / 2, width / 2, length / 2)));

  Common::AlwaysAccept<F> model;

  PET3D::Event<F> event(Point(0, 0, 0), Vector(1, 0, 0));
  RNG rng;

  SECTION("Two volumes") {

    typename PET3D::Full::Scanner<F, S>::FullResponse response;
    S hits = scanner.exact_detect(rng, model, event, response);
    CHECK(hits == 2);

    CHECK(response.detector1 == 0);
    CHECK(response.detector2 == 1);

    CHECK(response.d1_entry == VApprox(Point(R - height / 2, 0.0, 0.0)));
    CHECK(response.d1_deposition == VApprox(Point(R - height / 2, 0.0, 0.0)));
    CHECK(response.d1_exit == VApprox(Point(R + height / 2, 0.0, 0.0)));

    CHECK(response.d2_entry == VApprox(Point(-R + height / 2, 0.0, 0.0)));
    CHECK(response.d2_deposition == VApprox(Point(-R + height / 2, 0.0, 0.0)));
    CHECK(response.d2_exit == VApprox(Point(-R - height / 2, 0.0, 0.0)));

    CHECK(response.origin == VApprox(Point(0.0, 0.0, 0.0)));
  }
}

TEST_CASE("") {

  PET3D::Full::Scanner<F, S> scanner;

  F length = 0.500;
  F width = 0.007;
  F height = 0.019;
  F R = 0.500;

  Common::AlwaysAccept<F> model;

  PET3D::Event<F> event(Point(0, 0, 0), Vector(1, 0, 0));
  RNG rng;

  scanner.add_volume(PET3D::Full::BoxVolume<F>::AAB(
      Point(R / 2 - height / 2, -width / 2, -length / 2),
      Point(R / 2 + height / 2, width / 2, length / 2)));

  scanner.add_volume(PET3D::Full::BoxVolume<F>::AAB(
      Point(R - height / 2, -width / 2, -length / 2),
      Point(R + height / 2, width / 2, length / 2)));

  scanner.add_volume(PET3D::Full::BoxVolume<F>::AAB(
      Point(-R + height / 2, -width / 2, -length / 2),
      Point(-R - height / 2, width / 2, length / 2)));

  SECTION("Three volumes") {

    typename PET3D::Full::Scanner<F, S>::FullResponse response;
    S hits = scanner.exact_detect(rng, model, event, response);
    CHECK(hits == 2);

    CHECK(response.detector1 == 0);
    CHECK(response.detector2 == 2);

    CHECK(response.d1_entry == VApprox(Point(R / 2 - height / 2, 0.0, 0.0)));
    CHECK(response.d1_deposition ==
          VApprox(Point(R / 2 - height / 2, 0.0, 0.0)));
    CHECK(response.d1_exit == VApprox(Point(R / 2 + height / 2, 0.0, 0.0)));

    CHECK(response.d2_entry == VApprox(Point(-R + height / 2, 0.0, 0.0)));
    CHECK(response.d2_deposition == VApprox(Point(-R + height / 2, 0.0, 0.0)));
    CHECK(response.d2_exit == VApprox(Point(-R - height / 2, 0.0, 0.0)));

    CHECK(response.origin == VApprox(Point(0.0, 0.0, 0.0)));
  }
}

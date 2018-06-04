#include "util/test.h"

#include <random>

#include "generic_scanner.h"
#include "scanner_builder.h"
#include "common/types.h"
#include "common/model.h"

TEST("2d/barrel/generic_scanner/detection/") {

  using SquareDetector = PET2D::Barrel::SquareDetector<F>;
  using Detector = PET2D::Barrel::GenericScanner<SquareDetector, S, 32>;
  using Event = Detector::Event;
  using RNG = std::mt19937_64;
#if !_MSC_VER
  using FullResponse = typename Detector::FullResponse;
#else
  using FullResponse = Detector::FullResponse;
#endif

  Detector scanner = PET2D::Barrel::ScannerBuilder<Detector>::build_single_ring(
      0.43, 32, F(.005), F(.019));

  Common::AlwaysAccept<F> model;

  SECTION("central event") {
    Event e(0, 0, 0);
    FullResponse full_response;
    RNG rng;

    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 16);
    CHECK(full_response.lor.second == 0);
    CHECK(full_response.dl == Approx(0));
  }

  SECTION("non central event right") {
    Event e(0.2, 0.0, 0);
    FullResponse full_response;
    RNG rng;

    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 16);
    CHECK(full_response.lor.second == 0);
    CHECK(full_response.dl == Approx(0.43 + 0.2 - (0.43 - 0.2)));
  }

  SECTION("non central event left ") {
    Event e(-0.2, 0.0, 0);
    FullResponse full_response;
    RNG rng;

    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 16);
    CHECK(full_response.lor.second == 0);
    CHECK(full_response.dl == Approx(0.43 - 0.2 - (0.43 + 0.2)));
  }

  SECTION("non central event 90 degrees") {
    Event e(0.0, 0.2, M_PI / 2);
    FullResponse full_response;
    RNG rng;

    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 24);
    CHECK(full_response.lor.second == 8);
    CHECK(full_response.dl == Approx(0.43 + 0.2 - (0.43 - 0.2)));
  }

  SECTION("non central event 90 degrees oposite") {
    Event e(0.0, -0.2, M_PI / 2);
    FullResponse full_response;
    RNG rng;

    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 24);
    CHECK(full_response.lor.second == 8);
    CHECK(full_response.dl == Approx(0.43 - 0.2 - (0.43 + 0.2)));
  }
}

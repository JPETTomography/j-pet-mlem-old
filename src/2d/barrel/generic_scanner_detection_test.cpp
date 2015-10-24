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

  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  Detector scanner = PET2D::Barrel::ScannerBuilder<Detector>::build_single_ring(
      0.43, 32, F(.005), F(.019));

  SECTION("central event") {
    Event e(0, 0, 0);
    FullResponse full_response;
    RNG rng;
    Common::AlwaysAccept<F> model;
    auto hits = scanner.detect(rng, model, e, full_response);

    CHECK(hits == 2);

    CHECK(full_response.lor.first == 16);
    CHECK(full_response.lor.second == 0);
    CHECK(full_response.dl == Approx(0));
  }

}

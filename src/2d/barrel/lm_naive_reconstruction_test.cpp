
#include "util/test.h"

#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "common/types.h"
#include "common/model.h"
#include "2d/barrel/geometry.h"
#include "2d/barrel/lm_reconstruction.h"

TEST("2d/barrel/lm_reconstruction/naive") {

  using SquareDetector = PET2D::Barrel::SquareDetector<F>;
  using Detector = PET2D::Barrel::GenericScanner<SquareDetector, S, 32>;
  using Event = Detector::Event;
  using RNG = std::mt19937_64;

  using Response = typename Detector::Response;
  using FullResponse = typename Detector::FullResponse;

  Detector scanner = PET2D::Barrel::ScannerBuilder<Detector>::build_single_ring(
      0.43, 32, F(.005), F(.019));

  util::ibstream in_geometry("test_input/g_test");
  REQUIRE(in_geometry.is_open());

  PET2D::Barrel::Geometry<F, S> geometry(in_geometry);

  CHECK(geometry.n_detectors == 32);
  CHECK(geometry.grid.n_columns == 64);
  CHECK(geometry.grid.n_rows == 64);
  CHECK(geometry.grid.pixel_size == Approx(0.01));

  PET2D::Barrel::LMReconstruction<F, S,32> reconstruction(geometry, 0.04);

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

    Response response = scanner.response_wo_error(full_response);

    CHECK(response.lor.first == 16);
    CHECK(response.lor.second == 0);
    CHECK(response.dl == Approx(0));

    reconstruction.add(response);
  }
}

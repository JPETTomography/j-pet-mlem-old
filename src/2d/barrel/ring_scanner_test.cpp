#include <iostream>
#include <fstream>

#include "util/test.h"

#include "common/model.h"
#include "ring_scanner.h"
#include "scanner_builder.h"

#include "common/types.h"

using SquareDetector = PET2D::Barrel::SquareDetector<F>;
using Scanner = PET2D::Barrel::RingScanner<SquareDetector, S, 512>;
using Model = Common::AlwaysAccept<F>;

TEST("2d/barrel/scanner/math") {

  std::ifstream in("test_input/scanner_test.tab");

  if (!in) {
    WARN(
        "cannot open file `test_input/scanner_test.tab', "
        "evaluate `math/polygon_test.nb'");
    return;
  }

  double r, w, h;
  int n_detectors;
  int n_events;

  in >> r >> n_detectors >> w >> h >> n_events;

  Scanner ring = PET2D::Barrel::ScannerBuilder<Scanner>::build_single_ring(
      r, n_detectors, w, h);

  for (int i_event = 0; i_event < n_events; ++i_event) {
    double x, y, phi;
    in >> x >> y >> phi;

    Scanner::Event event(x, y, phi);

    in >> n_detectors;

    std::vector<int> detector(n_detectors);
    for (int i = 0; i < n_detectors; i++) {
      in >> detector[i];
      detector[i]--;  // mathematica counts positions from 1
    }

    for (int i = 0; i < n_detectors; i++) {
      in >> x >> y;
      Scanner::Point p1(x, y);

      in >> x >> y;
      Scanner::Point p2(x, y);

      auto inters = ring[detector[i]].intersections(event);
      CHECK(inters.size() == 2);

      CHECK(std::min(p1.x, p2.x) ==
            Approx(std::min(inters[0].x, inters[1].x)).epsilon(1e-4));
      CHECK(std::max(p1.x, p2.x) ==
            Approx(std::max(inters[0].x, inters[1].x)).epsilon(1e-4));
    }

    // this is not yet a complete tests....
    Model model;
    double position;
#if !_MSC_VER
    typename Scanner::Response response;
#else
    Scanner::Response response;
#endif
    auto hits = ring.detect(model, model, event, response);

    if (hits >= 2) {
      CHECK(std::find(detector.begin(), detector.end(), response.lor.first) !=
            detector.end());
      CHECK(std::find(detector.begin(), detector.end(), response.lor.second) !=
            detector.end());
    }
  }
}

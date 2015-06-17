#include <iostream>
#include <fstream>

#include "util/test.h"

#include "common/model.h"
#include "ring_scanner.h"
#include "scanner_builder.h"

using SquareDetector = PET2D::Barrel::SquareDetector<float>;
using Scanner = PET2D::Barrel::RingScanner<SquareDetector, 512, short>;
using Model = PET2D::Barrel::AlwaysAccept<float>;

TEST("2d/barrel/scanner/math") {

  std::ifstream in("math/scanner_test.tab");

  if (!in) {
    WARN("cannot open file `math/scanner_test.tab'");
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
            Approx(std::min(inters[0].x, inters[1].x)).epsilon(1e-5));
      CHECK(std::max(p1.x, p2.x) ==
            Approx(std::max(inters[0].x, inters[1].x)).epsilon(1e-5));
    }

    // this is not yet a complete tests....
    Model model;
    double position;
    typename Scanner::Response response;
    auto hits = ring.detect(model, model, event, response);

    if (hits >= 2) {
      CHECK(std::find(detector.begin(), detector.end(), response.lor.first) !=
            detector.end());
      CHECK(std::find(detector.begin(), detector.end(), response.lor.second) !=
            detector.end());
    }
  }
}

#include <iostream>
#include <fstream>

#include "util/test.h"

#include "model.h"
#include "detector_ring.h"

using namespace PET2D;
using namespace PET2D::Barrel;

TEST("2d/barrel/detector_ring/math") {

  std::ifstream in("math/detector_ring_test.tab");

  if (!in) {
    WARN("cannot open file `math/detector_ring_test.tab'");
    return;
  }

  double r, w, h;
  int n_detectors;
  int n_events;

  in >> r >> n_detectors >> w >> h >> n_events;

  DetectorRing<> ring(n_detectors, r, w, h);

  for (int i_event = 0; i_event < n_events; ++i_event) {
    double x, y, phi;
    in >> x >> y >> phi;

    DetectorRing<>::Event event(x, y, phi);

    in >> n_detectors;

    std::vector<int> detector(n_detectors);
    for (int i = 0; i < n_detectors; i++) {
      in >> detector[i];
      detector[i]--;  // mathematica counts positions from 1
    }

    for (int i = 0; i < n_detectors; i++) {
      double x, y;

      in >> x >> y;
      DetectorRing<>::Point p1(x, y);

      in >> x >> y;
      DetectorRing<>::Point p2(x, y);

      auto inters = ring[detector[i]].intersections(event);
      CHECK(inters.size() == 2);

      CHECK(std::min(p1.x, p2.x) ==
            Approx(std::min(inters[0].x, inters[1].x)).epsilon(1e-13));
      CHECK(std::max(p1.x, p2.x) ==
            Approx(std::max(inters[0].x, inters[1].x)).epsilon(1e-13));
    }

    // this is not yet a complete tests....
    DetectorRing<>::LOR lor;
    AlwaysAccept<> model;
    double position;
    auto hits = ring.detect(model, model, event, lor, position);

    if (hits >= 2) {
      CHECK(std::find(detector.begin(), detector.end(), lor.first) !=
            detector.end());
      CHECK(std::find(detector.begin(), detector.end(), lor.second) !=
            detector.end());
    }
  }
}

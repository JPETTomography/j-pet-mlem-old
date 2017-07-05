#ifndef GATE_SCANNER_BUILDER_H
#define GATE_SCANNER_BUILDER_H

#include "2d/gate/gate_volume.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/square_detector.h"

namespace Gate {
namespace D2 {

template <typename FType, typename SType> class GenericScannerBuilder {
 public:
  using F = FType;
  using S = SType;
  using Volume = Gate::D2::Volume<F>;
  using Detector = PET2D::Barrel::SquareDetector<F>;
  using Scanner = PET2D::Barrel::GenericScanner<Detector, S>;
  void build(Volume* vol, Scanner* scanner) {

    if (vol->is_sd()) {
      auto box = dynamic_cast<Gate::D2::Box<F>*>(vol);
      if (box) {
        Detector d(box->lengthX, box->lengthY, 1);
        d += box->translation();
        scanner->push_back(d);
      } else {
        fprintf(stderr, "unsupported volume");
      }
    }
    for (auto daughter = vol->daughters(); daughter != vol->daughters_end();
         daughter++) {
      build(*daughter, scanner);
    }
  }
};
}
}

#endif  // GATE_SCANNER_BUILDER_H

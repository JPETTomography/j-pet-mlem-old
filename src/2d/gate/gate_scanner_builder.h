#ifndef GATE_SCANNER_BUILDER_H
#define GATE_SCANNER_BUILDER_H

#include "2d/gate/gate_volume.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/square_detector.h"

namespace Gate {
namespace D2 {

template <typename FType, typename SType> class GenericScannerBuilder {
  using F = FType;
  using S = SType;
  using Volume = Gate::D2::Volume<F>;
  using Scanner =
      PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S>;
  void build(Volume* vol, Scanner& scanner){};
};
}
}

#endif  // GATE_SCANNER_BUILDER_H

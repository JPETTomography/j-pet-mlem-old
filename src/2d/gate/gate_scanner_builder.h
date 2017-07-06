#ifndef GATE_SCANNER_BUILDER_H
#define GATE_SCANNER_BUILDER_H

#include "2d/gate/gate_volume.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/square_detector.h"
#include "2d/geometry/transformation.h"

namespace Gate {
namespace D2 {

template <typename FType, typename SType> class GenericScannerBuilder {
 public:
  using F = FType;
  using S = SType;
  using Volume = Gate::D2::Volume<F>;
  using Detector = PET2D::Barrel::SquareDetector<F>;
  using Scanner = PET2D::Barrel::GenericScanner<Detector, S>;
  using Vector = typename Detector::Vector;
  using Transformation = PET2D::Transformation<F>;
  using Repeater = Gate::D2::Repeater<F>;

  void build(Volume* vol, Scanner* scanner, Transformation transformation) {
    // Transformation local_transform(vol->rotation(), vol->translation());
    Transformation local_transform = vol->transformation();
    Repeater* repeater;
    if ((repeater = vol->repeater()) != nullptr) {
      vol->detach_repeater();
      for (int i = 0; i < repeater->n; i++) {
        vol->set_transformation(
            new Transformation(local_transform * repeater->operator[](i)));
        build(vol, scanner, transformation);
      }
    } else {

      if (vol->is_sd()) {
        auto box = dynamic_cast<Gate::D2::Box<F>*>(vol);
        if (box) {
          Detector d(box->lengthX, box->lengthY, 1);
          d.transform(local_transform);
          d.transform(transformation);
          scanner->push_back(d);
        } else {
          fprintf(stderr, "unsupported volume");
        }
      }
    }
    for (auto daughter = vol->daughters(); daughter != vol->daughters_end();
         daughter++) {
      build(*daughter, scanner, transformation * local_transform);
    }
  }

  void build(Volume* vol, Scanner* scanner) {
    build(vol, scanner, Transformation());
  }
};
}
}

#endif  // GATE_SCANNER_BUILDER_H

#ifndef GATE_SCANNER_BUILDER_H
#define GATE_SCANNER_BUILDER_H

#include "2d/gate/gate_volume.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/square_detector.h"
#include "2d/geometry/transformation.h"

namespace Gate {
namespace D2 {

template <typename FType, typename SType, int MaxDetectors = MAX_DETECTORS>
class GenericScannerBuilder {
 public:
  using F = FType;
  using S = SType;
  using Volume = Gate::D2::Volume<F>;
  using Detector = PET2D::Barrel::SquareDetector<F>;
  using Scanner = PET2D::Barrel::GenericScanner<Detector, S, MaxDetectors>;
  using Vector = typename Detector::Vector;
  using Transformation = PET2D::Transformation<F>;
  using Repeater = Gate::D2::Repeater<F>;

  void build(const Volume* vol,
             Scanner* scanner,
             Transformation transformation) {
    Transformation local_transform = vol->transformation();
    Repeater* repeater;
    if ((repeater = vol->repeater()) != nullptr) {
      auto repeater = const_cast<Volume*>(vol)->detach_repeater();
      for (int i = 0; i < repeater->n; i++) {
        const_cast<Volume*>(vol)->set_transformation(
            new Transformation(repeater->operator[](i)*local_transform));
        build(vol, scanner, transformation);
      }
      const_cast<Volume*>(vol)->attach_repeater(repeater);
      const_cast<Volume*>(vol)
          ->set_transformation(new Transformation(local_transform));
    } else {

      if (vol->is_sd()) {
        auto box = dynamic_cast<const Gate::D2::Box<F>*>(vol);
        if (box) {
          Detector d(box->lengthX, box->lengthY, 1);
          d.transform(local_transform);
          d.transform(transformation);
          scanner->push_back(d);
        } else {
          fprintf(stderr, "unsupported volume");
        }
      }

      for (auto daughter = vol->daughters(); daughter != vol->daughters_end();
           daughter++) {
        build(*daughter, scanner, transformation * local_transform);
      }
    }
  }

  void build(Volume* vol, Scanner* scanner) {
    build(vol, scanner, Transformation());
  }

  S count_cristals(const Volume* vol, S counter) {
    S local_counter = counter;
    Repeater* repeater;
    if ((repeater = vol->repeater()) != nullptr) {
      auto repeater = const_cast<Volume*>(vol)->detach_repeater();
      for (int i = 0; i < repeater->n; i++)
        local_counter = count_cristals(vol, local_counter);
      const_cast<Volume*>(vol)->attach_repeater(repeater);
    } else {

      if (vol->is_sd()) {
        local_counter++;
      }

      for (auto daughter = vol->daughters(); daughter != vol->daughters_end();
           daughter++)
        local_counter = count_cristals(*daughter, local_counter);
    }
    return local_counter;
  }

  S count_cristals(const Volume* vol) { return count_cristals(vol, 0); }

  Scanner build_with_8_symmetries(Volume* vol){

  };
};
}
}

#endif  // GATE_SCANNER_BUILDER_H

#include "util/test.h"

#include "2d/gate/gate_volume.h"
#include "2d/gate/gate_scanner_builder.h"
#include "common/types.h"

TEST("2d Gate volume") {
  SECTION("Constructing") {
    Gate::D2::Box<float> world;
    Gate::D2::GenericScannerBuilder<F, S> builder;
  }
}

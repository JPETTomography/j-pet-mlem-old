#include <fstream>

#include "util/test.h"

#include "detector_set.h"

#include "scanner_builder.h"
#include "generic_scanner.h"
#include "square_detector.h"

#include "common/types.h"

using Builder = PET2D::Barrel::ScannerBuilder<
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S, 128>>;

TEST("Serialise") {
  auto scanner = Builder::build_single_ring(200.0, 8, 0.007, 0.019);

  std::ofstream out("det_set_8.json");

  scanner.serialize(out);
  out.close();

  std::ifstream in("det_set_8.json");
  auto scanner_copy = Builder::deserialize(in);

  for (int d = 0; d < scanner.size(); d++) {
    REQUIRE(scanner[d].approx_equal_dihedral(scanner_copy[d], 1e-5));
  }
  in.close();
}

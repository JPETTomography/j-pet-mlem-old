#include <fstream>

#include "util/test.h"

#include "detector_set.h"

#include "scanner_builder.h"
#include "generic_scanner.h"
#include "square_detector.h"

#include "common/types.h"

using Scanner =
    PET2D::Barrel::GenericScanner<PET2D::Barrel::SquareDetector<F>, S, 128>;
using Builder = PET2D::Barrel::ScannerBuilder<Scanner>;

TEST("Serialise") {
  F R = 200.0;
  F w = 0.007;
  F h = 0.019;
  auto scanner = Builder::build_single_ring(R, 8, w, h);

  SECTION("serialize dets") {
    std::string name("det_set_8.txt");
    std::ofstream out(name);
    scanner.serialize_detectors(out);
    out.close();

    std::ifstream in(name);
    auto scanner_copy = Builder::deserialize(in);

    for (int d = 0; d < scanner.size(); d++) {
      REQUIRE(scanner[d].approx_equal_dihedral(scanner_copy[d], 1e-5));
    }
    in.close();
    std::remove(name.c_str());
  }

  SECTION("serialize") {
    std::string name_d("det_set_8.txt");
    std::ofstream out_d(name_d);
    std::string name_s("det_set_8_sym.txt");
    std::ofstream out_s(name_s);
    scanner.serialize(out_d, out_s);
    out_d.close();
    out_s.close();

    auto descriptor = scanner.symmetry_descriptor();

    std::ifstream in_d(name_d);
    auto scanner_copy = Builder::deserialize(in_d);
    in_d.close();
    std::remove(name_d.c_str());
    std::ifstream in_s(name_s);
    scanner_copy.set_symmetry_descriptor(
        PET2D::Barrel::SymmetryDescriptor<S>::deserialize(in_s));
    in_s.close();
    std::remove(name_s.c_str());
    for (int d = 0; d < scanner.size(); d++) {
      REQUIRE(scanner[d].approx_equal_dihedral(scanner_copy[d], 1e-5));
    }
    auto descriptor_copy = scanner_copy.symmetry_descriptor();
    REQUIRE(descriptor.n_detectors == descriptor_copy.n_detectors);
    REQUIRE(descriptor.n_symmetries == descriptor_copy.n_symmetries);
    for (S d = 0; d < descriptor.n_detectors; d++) {
      for (S s = 0; s < descriptor.n_symmetries; s++) {
        REQUIRE(descriptor.symmetric_detector(d, s) ==
                descriptor_copy.symmetric_detector(d, s));
      }
    }
  }

  SECTION("radius") {
    auto rs = scanner.min_max_radius();
    auto r_min = std::sqrt(R * R + w * w / 4);
    auto r_max = std::sqrt((R + h) * (R + h) + w * w / 4);
    REQUIRE(rs.first == Approx(r_min));
    REQUIRE(rs.second == Approx(r_max));
  }
}

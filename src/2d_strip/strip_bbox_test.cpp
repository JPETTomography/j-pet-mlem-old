#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "reconstruction.h"

StripDetector<double> detector(500, 1000, 200, 200, 5, 5, 10, 63);
Reconstruction<double> reconstructor(1, detector);

void check(double angle, double bby_value, double bbz_value) {

  double inv_pow_sigma_dl =
      (1.0 / (detector.get_sigma_dl() * detector.get_sigma_dl()));
  double inv_pow_sigma_z =
      (1.0 / (detector.get_sigma_z() * detector.get_sigma_z()));

  double A = ((4.0 / (cos(angle) * cos(angle))) * inv_pow_sigma_dl) +
             (2.0 * tan(angle) * tan(angle) * inv_pow_sigma_z);
  double B = -4.0 * tan(angle) * inv_pow_sigma_z;
  double C = 2.0 * inv_pow_sigma_z;
  double B_2 = (B / 2.0) * (B / 2.0);

  CHECK(reconstructor.bby(A, C, B_2) == Approx(bby_value).epsilon(1e-4));
  CHECK(reconstructor.bbz(A, C, B_2) == Approx(bbz_value).epsilon(1e-4));
}

TEST_CASE("bounding box test", "4") {

  check(0.0470448, 94.3954, 21.6737);
  check(-0.594145, 78.3053, 56.9959);
  check(0.20029, 92.6108, 28.3458);
  check(-0.571667, 79.4745, 55.3539);
  check(-0.420542, 86.266, 44.0276);
}

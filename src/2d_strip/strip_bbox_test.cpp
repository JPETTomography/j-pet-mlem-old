#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include "strip_detector.h"

TEST_CASE("strip/bbox") {

  StripDetector<> detector(500, 1000, 200, 200, 5, 5, 10, 63);

  struct {
    double angle;
    double bby_value;
    double bbz_value;
  } v[] = { { 0.0470448, 94.3954, 21.6737 },
            { -0.594145, 78.3053, 56.9959 },
            { 0.20029, 92.6108, 28.3458 },
            { -0.571667, 79.4745, 55.3539 },
            { -0.420542, 86.266, 44.0276 } };

  for (int i = 0; i < sizeof(v) / sizeof(*v); ++i) {
    auto inv_pow_sigma_dl = (1.0 / (detector.sigma_dl * detector.sigma_dl));
    auto inv_pow_sigma_z = (1.0 / (detector.sigma_z * detector.sigma_z));
    auto angle = v[i].angle;
    auto bby_value = v[i].bby_value;
    auto bbz_value = v[i].bbz_value;

    auto A = ((4.0 / (cos(angle) * cos(angle))) * inv_pow_sigma_dl) +
             (2.0 * tan(angle) * tan(angle) * inv_pow_sigma_z);
    auto B = -4.0 * tan(angle) * inv_pow_sigma_z;
    auto C = 2.0 * inv_pow_sigma_z;
    auto B_2 = (B / 2.0) * (B / 2.0);

    CHECK(detector.bb_y(A, C, B_2) == Approx(bby_value).epsilon(1e-4));
    CHECK(detector.bb_z(A, C, B_2) == Approx(bbz_value).epsilon(1e-4));
  }
}

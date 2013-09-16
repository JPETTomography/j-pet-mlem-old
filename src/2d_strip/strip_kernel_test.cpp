#include "catch.hpp"
#include <vector>
#include <cmath>

#include "util/png_writer.h"
#include "util/bstream.h"
#include "util/svg_ostream.h"

#include"reconstruction.h"

TEST_CASE("kernel tests", "kernel") {


Reconstruction<> reconstructor(1,500,1000,200,5,10,63);
Reconstruction<>::Point  from_center(0,0);
double y=0.0;
double tan_value=0.0;
double inv_cos=1.0;
double pow_inv_cos=1.0;

double value=reconstructor.kernel(y, tan_value, inv_cos, pow_inv_cos, from_center);
CHECK(value == Approx(1.1372205719261035e-7).epsilon(1e-14));
}

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "common/types.h"

#include "2d/strip/kernel_monte_carlo.h"
#include "2d/geometry/pixel_map.h"
#include "2d/geometry/pixel.h"

#include "2d/strip/options.h"
#include "2d/strip/response.h"
#include "2d/strip/reconstruction.h"
#include "2d/strip/gaussian_kernel.h"

using Pixel = PET2D::Pixel<int>;
using Output = PET2D::PixelMap<Pixel, F>;

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Reconstruction = PET2D::Strip::Reconstruction<F, Kernel>;
using Scanner = PET2D::Strip::Scanner<F, S>;

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  cl.parse_check(argc, argv);

  strip_integral();
  strip_integral_theta();
  strip_integral_event();
  strip_integral_theta_event();
  strip_gauss_kernel();
  strip_gauss_kernel_integral();
}

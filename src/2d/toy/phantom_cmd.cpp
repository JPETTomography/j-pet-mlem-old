
#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/backtrace.h"

#include "2d/geometry/phantom.h"
#include "common/phantom_monte_carlo.h"
#include "2d/toy/gauss_scanner.h"
#include "common/types.h"

using RNG = util::random::tausworthe;
using Scanner = PET2D::Toy::GaussScanner<F>;
using Phantom = PET2D::Phantom<RNG, F>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;

main(int argc, const char* argv[]) {

  CMDLINE_TRY
  cmdline::parser cl;

  cl.add<double>("length-x", 'l', "length in x", false, 1);
  cl.add<double>("length-y", '\0', "length in y", false, 1);

  cl.add<double>("s-pixel", 'p', "pixel size", false, 0.01);
  cl.add<int>("n-pixels", 'n', "number of pixels", cmdline::dontsave, 0);
  cl.add<int>("n-z-pixels", 0, "number of z pixels", false);
  cl.add<int>("n-y-pixels", 0, "number of y pixels", false);
  cl.add<double>(
      "s-z", 0, "TOF sigma along z axis", cmdline::alwayssave, 0.015);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
  cl.add<int>("emissions", 'e', "number of emissions", false, 0);
  cl.add<double>("scale", '\0', "scale factor", false, 1);

  Phantom phantom(cl.get<double>("scale"));
  for (auto& fn : cl.rest()) {
    std::ifstream in_phantom(fn);
    phantom << in_phantom;
  }

  phantom.calculate_cdf();
  Scanner scanner(cl.get<double>("s-z"), cl.get<double>("s-dl"));
  MonteCarlo monte_carlo(phantom, scanner);

  CMDLINE_CATCH
}

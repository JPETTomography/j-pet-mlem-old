#include "cmdline.h"
#include "util/backtrace.h"
#include "util/cmdline_hooks.h"
#include "util/cmdline_types.h"

#include "2d/gate/gate_scanner_builder.h"
#include "2d/gate/gate_volume_builder.h"

#include "common/types.h"

int main(int argc, char* argv[]) {
  using Volume = Gate::D2::Volume<F>;

  CMDLINE_TRY

  cmdline::parser cl;

  cl.add<cmdline::path>("output", 'o', "ouput file", true);
  cl.add("big-barrel", 'B', "creates big barrel description");
  cl.add("small-barrel", 'b', "creates small barrel description");
  cl.add("new-module", 'n', "new modules description");
  cl.add("full", 'f', "full detector: new + big barrel");
  cl.add("ideal", 'i', "ideal geometry with 1 layer and 384 strips ");

  cl.parse_check(argc, argv);

  auto out = cl.get<cmdline::path>("output");
  auto extension = out.ext();
  if (extension.empty())
    extension = "txt";
  auto wo_ext = out.wo_ext();

  std::ofstream o_dets(cl.get<cmdline::path>("output") + "_dets." + extension);
  std::ofstream o_syms(cl.get<cmdline::path>("output") + "_syms." + extension);

  Volume* world;
  if (cl.exist("big-barrel")) {
    world = Gate::D2::build_big_barrel_volume<F>();
  } else if (cl.exist("full")) {
    world = Gate::D2::build_new_full_scanner_volume<F>();
  } else if (cl.exist("ideal")) {
    world = Gate::D2::build_ideal_scanner_volume<F>();
  }

  const int n_detectors = Gate::D2::count_cristals<F, S>(world);
  if (n_detectors <= 256) {
    auto builder = new Gate::D2::GenericScannerBuilder<F, S, 256>;
    auto scanner = builder->build_with_8_symmetries(world);

    scanner.serialize(o_dets, o_syms);
    o_dets.close();
    o_syms.close();
    return 0;
  } else if (n_detectors <= 512) {
    auto builder = new Gate::D2::GenericScannerBuilder<F, S, 512>;
    auto scanner = builder->build_with_8_symmetries(world);

    scanner.serialize(o_dets, o_syms);

    return 0;
  }

  o_dets.close();
  o_syms.close();

  CMDLINE_CATCH
  return 0;
}

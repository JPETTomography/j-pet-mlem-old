#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/backtrace.h"

int main(int argc, char* argv[]) {

  CMDLINE_TRY

  cmdline::parser cl;

  cl.add<std::string>("output", 'o', "ouput file", true);
  cl.add("--big-barrel", 'B', "creates big barrel description");
  cl.add("--small-barrel", 'b', "creates small barrel description");
  cl.add("--new-module", 'n', "new modules description");
  cl.add("--full", 'f', "full detector: new + big barrel");

  CMDLINE_CATCH
  return 0;
}

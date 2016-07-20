#include "test.h"
#include "cmdline.h"
#include "cmdline_types.h"

TEST("util/cmdline") {
  SECTION("filename") {
    cmdline::path fn("/tmp/dir/dir.ext/file.name_some.ext");
    CHECK(fn.ext() == ".ext");
    CHECK(fn.wo_path() == "file.name_some.ext");
  }
  SECTION("index") {
    cmdline::path fn("/tmp/dir/dir.ext/file.name_some.ext_001");
    CHECK(fn.scan_index() == 1);
  }
}

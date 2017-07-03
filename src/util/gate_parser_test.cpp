#include "test.h"
#include "util/gate_parser.h"

TEST("util/gate_parser") {
  SECTION("one arg parse") {
    Gate::Parser parser;
    auto command_chain = parser.parse("/gate/world/name box");
    auto command = command_chain.commands();

    CHECK("gate" == *command);
    command++;
    CHECK("world" == *command);
    command++;
    CHECK("name" == *command);
  }
}

#include "test.h"
#include "util/gate_parser.h"

TEST("util/gate_parser") {
  SECTION("zero arg parse") {
    Gate::Parser parser;
    auto command_chain = parser.parse(" /gate/world/name ");

    CHECK(command_chain.command_length() == 3);
    CHECK(command_chain.n_arguments() == 0);

    auto command = command_chain.commands();

    CHECK("gate" == *command);
    command++;
    CHECK("world" == *command);
    command++;
    CHECK("name" == *command);
  }

  SECTION("one arg parse") {
    Gate::Parser parser;
    auto command_chain = parser.parse("/gate/world/name box");
    auto command = command_chain.commands();

    CHECK("gate" == *command);
    command++;
    CHECK("world" == *command);
    command++;
    CHECK("name" == *command);

    auto argument = command_chain.arguments();

    CHECK("box" == *argument);
  }

  SECTION("four arg parse") {
    Gate::Parser parser;
    auto command_chain =
        parser.parse("/gate/box/setTranslation 0.0 1.0 2.0 cm");
    auto command = command_chain.commands();

    CHECK(command_chain.command_length() == 3);
    CHECK(command_chain.n_arguments() == 4);

    CHECK("gate" == *command);
    command++;
    CHECK("box" == *command);
    command++;
    CHECK("setTranslation" == *command);

    auto argument = command_chain.arguments();

    CHECK("0.0" == *argument);
    argument++;
    CHECK("1.0" == *argument);
    argument++;
    CHECK("2.0" == *argument);
    argument++;
    CHECK("cm" == *argument);
  }
}

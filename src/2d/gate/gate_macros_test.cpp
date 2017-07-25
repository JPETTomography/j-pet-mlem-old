#include <iostream>
#include <fstream>
#include <sstream>

#include "util/test.h"
#include "gate_volume.h"
#include "gate_macros.h"
#include "common/types.h"

TEST("Gate macros", "only world") {
  using namespace Gate::D2;

  GateMacrosInterpreter<F, S> interpreter;
  std::string macro_name("src/2d/gate/world.txt");

  std::string macro_content(
      "/gate/geometry/setMaterialDatabase   ./GateMaterials.db\n");
  macro_content += +"/gate/world/geometry/setXLength 2. m\n";
  macro_content += "/gate/world/geometry/setYLength 2. m\n";
  macro_content += "/gate/world/geometry/setZLength 1. m";
  std::stringstream macro(macro_content);

  if (macro) {
    interpreter.interpret(macro);
    auto vol = interpreter.volume();
    auto world = dynamic_cast<Box<F>*>(vol);

    REQUIRE(world->lengthX == Approx(2.0));
    REQUIRE(world->lengthY == Approx(1.0));

  } else {
    WARN("Cannot open string stream");
  }
}

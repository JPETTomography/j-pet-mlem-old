#ifndef GATE_MACROS_H
#define GATE_MACROS_H

#include "gate_volume.h"
#include "util/gate_parser.h"

namespace Gate {
namespace D2 {

template <typename FType, typename SType> class GateMacrosInterpreter {
 public:
  using F = FType;
  void interpret(const std::istream& macro){};
  Gate::D2::Volume<F>* volume() { return 0; };
};
}
}

#endif  // GATE_MACROS_H

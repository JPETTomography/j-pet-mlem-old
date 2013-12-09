#include <fstream>
#include <iostream>

#include "cmdline_hooks.h"

namespace cmdline {

int loading = 0;

bool load(cmdline::parser& parser,
          string& value __attribute__((unused)),
          const std::string& arg) {
  std::ifstream in(arg);
  if (parser.exist("verbose")) {
    std::cout << "load: " << arg << std::endl;
  }
  if (!in.is_open()) {
    throw("cannot open input config file: " + arg);
  }
  ++loading;
  parser.try_parse(in, false, arg.c_str());
  --loading;
  return true;
}

bool not_from_file(cmdline::parser& parser __attribute__((unused)),
                   int& value __attribute__((unused)),
                   const std::string& arg __attribute__((unused))) {
  return loading == 0;
}

}  // cmdline

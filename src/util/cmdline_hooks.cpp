#include <fstream>
#include <iostream>

#include "cmdline_hooks.h"

namespace cmdline {

std::vector<std::string> dir_stack;

bool load(cmdline::parser& parser,
          string& value __attribute__((unused)),
          const std::string& arg) {
  std::string path = arg;
  if (!dir_stack.empty() && path.substr(0, 1) != "/") {
    path = dir_stack.back() + path;
  }
  if (parser.exist("verbose")) {
    std::cout << "load: " << path << std::endl;
  }
  std::ifstream in(path);
  if (!in.is_open()) {
    throw("cannot open input config file: " + path);
  }
  auto dir_sep = path.find_last_of('/');
  if (dir_sep != std::string::npos) {
    dir_stack.push_back(path.substr(0, dir_sep + 1));
  } else {
    dir_stack.push_back(std::string());
  }
  parser.try_parse(in, false, path.c_str());
  dir_stack.pop_back();
  return true;
}

bool not_from_file(cmdline::parser& parser __attribute__((unused)),
                   int& value __attribute__((unused)),
                   const std::string& arg __attribute__((unused))) {
  return dir_stack.empty();
}

}  // cmdline

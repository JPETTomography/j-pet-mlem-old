#include <fstream>
#include <iostream>

#include "cmdline_hooks.h"

namespace cmdline {

  static std::vector<std::string> dir_stack;

  bool load(cmdline::parser& parser, path& value, const std::string& arg) {
    (void)value;  // unused
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

  template <typename T>
  bool not_from_file(cmdline::parser& parser,
                     T& value,
                     const std::string& arg) {
    (void)parser, (void)value, (void)arg;  // unused
    return dir_stack.empty();
  }

  template bool not_from_file(cmdline::parser&, int&, std::string const&);
  template bool not_from_file(cmdline::parser&,
                              cmdline::path&,
                              std::string const&);

  void load_accompanying_config(cmdline::parser& parser, bool only_one) {
    // load config files accompanying phantom files
    for (cmdline::path fn : parser.rest()) {
      std::ifstream in(fn.wo_ext() + ".cfg");
      // check if file exists
      if (in.is_open()) {
        in.close();
        cmdline::load(parser, fn, fn.wo_ext() + ".cfg");
        if (only_one) {
          break;
        }
      }
    }
  }

}  // cmdline

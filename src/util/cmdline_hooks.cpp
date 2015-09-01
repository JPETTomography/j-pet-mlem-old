#include <fstream>
#include <iostream>

#include <unistd.h>
#include <sys/syslimits.h>

#include "cmdline_hooks.h"

namespace cmdline {

namespace {
static std::vector<std::string> dir_stack;
static cmdline::path exe_dir;
}

#if defined(_WIN32) || defined(WIN32)
#define DIR_SEP '\\'
#else
#define DIR_SEP '/'
#endif

void init(const char* argv0) {
  cmdline::path cmd(argv0);
  char buf[PATH_MAX];
  if (!getcwd(buf, PATH_MAX))
    return;
  cmdline::path cwd(buf);
  // remove relative components when starts with ../
  while (cmd.length() > 3 && cmd[0] == '.' && cmd[1] == '.' &&
         cmd[2] == DIR_SEP) {
    auto cwd_dir_sep = cwd.find_last_of(DIR_SEP);
    if (cwd_dir_sep == std::string::npos)
      return;
    cwd = cwd.substr(0, cwd_dir_sep);
    cmd = cmd.substr(3);
  }
  auto cmd_dir_sep = cmd.find_last_of(DIR_SEP);
  if (cmd_dir_sep != std::string::npos) {
    exe_dir = cwd + DIR_SEP + cmd.substr(0, cmd_dir_sep) + DIR_SEP;
  } else {
    exe_dir = cwd + DIR_SEP;
  }
}

bool load(cmdline::parser& parser, path& value, const std::string& arg) {
  (void)value;  // unused
  std::string path = arg;
  // if path is relative, use dir stack
  if (!dir_stack.empty() && path[0] != DIR_SEP) {
    path = dir_stack.back() + path;
  }
  auto verbose = parser.exist("verbose");
  if (verbose) {
    std::cout << "# load: " << path << std::endl;
  }
  std::ifstream in(path);
  if (!in.is_open()) {
    // if we cannot load it, try to load it relatively to exe dir config/arg.cfg
    if (exe_dir.length() && arg[0] != DIR_SEP) {
      auto fallback = exe_dir + "config" + DIR_SEP + arg + ".cfg";
      in.open(fallback);
      if (!in.is_open()) {
        throw("cannot open input config file: " + path);
      } else {
        path = fallback;
        if (verbose) {
          std::cout << "# falling back to: " << path << std::endl;
        }
      }
    } else {
      throw("cannot open input config file: " + path);
    }
  }
  auto dir_sep = path.find_last_of(DIR_SEP);
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
bool not_from_file(cmdline::parser& parser, T& value, const std::string& arg) {
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

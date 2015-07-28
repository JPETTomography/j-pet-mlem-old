#pragma once

#include "cmdline.h"
#include "cmdline_types.h"

/// Extensions for cmdline library

namespace cmdline {

/// Loads serialized command line parameters from config file given as
/// argument
bool load(cmdline::parser& parser, path& value, const std::string& arg);

/// Skips setting option if it comes from config file
template <typename T>
bool not_from_file(cmdline::parser& parser, T& value, const std::string& arg);

/// Loads accompanying `.cfg` files from paths given as command line arguments
void load_accompanying_config(cmdline::parser& parser, bool only_one = false);

}  // cmdline

#define CMDLINE_TRY try {

#define CMDLINE_CATCH                                                 \
  return 0;                                                           \
  }                                                                   \
  catch (cmdline::exception & ex) {                                   \
    if (ex.help()) {                                                  \
      std::cerr << ex.usage();                                        \
    }                                                                 \
    for (auto& msg : ex.errors()) {                                   \
      auto name = ex.name();                                          \
      if (name) {                                                     \
        std::cerr << "error at " << name << ": " << msg << std::endl; \
      } else {                                                        \
        std::cerr << "error: " << msg << std::endl;                   \
      }                                                               \
    }                                                                 \
  }                                                                   \
  catch (const std::string& ex) {                                     \
    std::cerr << "error: " << ex << std::endl;                        \
    util::print_backtrace(std::cerr);                                 \
  }                                                                   \
  catch (const char* ex) {                                            \
    std::cerr << "error: " << ex << std::endl;                        \
    util::print_backtrace(std::cerr);                                 \
  }                                                                   \
  catch (const std::exception& ex) {                                  \
    std::cerr << "error: " << ex.what() << std::endl;                 \
    util::print_backtrace(std::cerr);                                 \
  }                                                                   \
  return 1;

#pragma once

#include "cmdline.h"
#include "cmdline_types.h"

namespace cmdline {

/// loads serialized command line parameters from config file given as
/// argument
bool load(cmdline::parser& parser, path& value, const std::string& arg);

/// skips setting option if it comes from config file
template <typename T>
bool not_from_file(cmdline::parser& parser, T& value, const std::string& arg);

/// loads accompanying .cfg files from paths given as command line arguments
void load_accompanying_config(cmdline::parser& parser, bool only_one = false);

}  // cmdline

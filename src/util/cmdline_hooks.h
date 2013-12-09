#pragma once

#include "cmdline.h"
#include "cmdline_types.h"

namespace cmdline {

bool load(cmdline::parser& parser, string& value, const std::string& arg);
bool not_from_file(cmdline::parser& parser, int& value, const std::string& arg);

}  // cmdline

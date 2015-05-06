/// \page cmd_2d_strip 2D Strip
/// \brief 2D Strip Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_2d_strip_phantom
/// - \subpage cmd_2d_strip_reconstruction

#pragma once

#include "cmdline.h"

namespace PET2D {
namespace Strip {

/// Adds scanner specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_strip_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_strip_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// calculates all empty values from existing other parameters
void calculate_scanner_options(cmdline::parser& parser);

/// provides initialization list for creating scanner
#define __PET2D_STRIP(...) __VA_ARGS__  // just pass-through
#define PET2D_STRIP_SCANNER_CL(cl)            \
  __PET2D_STRIP(cl.get<double>("r-distance"), \
                cl.get<double>("s-length"),   \
                cl.get<int>("n-y-pixels"),    \
                cl.get<int>("n-z-pixels"),    \
                cl.get<double>("p-size"),     \
                cl.get<double>("p-size"),     \
                cl.get<double>("s-z"),        \
                cl.get<double>("s-dl"))

}  // Strip
}  // PET2D

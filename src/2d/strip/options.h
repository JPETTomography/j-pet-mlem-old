/// \page cmd_2d_strip 2D Strip
/// \brief 2D Strip Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_2d_strip_phantom
/// - \subpage cmd_2d_strip_reconstruction
///
/// Workflow
/// --------
/// \f[
///  \\ \mathit{phantom~desc.} \rightarrow \mathtt{2d\_strip\_phantom}
///                            \rightarrow \mathit{response}
///                            \rightarrow \mathtt{2d\_strip\_reconstruction}
/// \f]
///
/// \note See \ref cmd_2d_strip_phantom and \ref cmd_2d_strip_reconstruction for
/// usage examples.

#pragma once

#include "cmdline.h"
#include "common/options.h"

namespace PET2D {
namespace Strip {

/// Adds scanner specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_strip_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_strip_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// calculates all empty values from existing other parameters
void calculate_scanner_options(cmdline::parser& parser, int argc = 1);

/// provides initialization list for creating scanner
#define __PET2D_STRIP(...) __VA_ARGS__  // just pass-through
#define PET2D_STRIP_SCANNER_CL(cl)                        \
  __PET2D_STRIP(cl.get<std::vector<double>>("radius")[0], \
                cl.get<double>("length"),                 \
                cl.get<int>("n-y-pixels"),                \
                cl.get<int>("n-z-pixels"),                \
                cl.get<double>("s-pixel"),                \
                cl.get<double>("s-pixel"),                \
                cl.get<double>("s-z"),                    \
                cl.get<double>("s-dl"))

}  // Strip
}  // PET2D

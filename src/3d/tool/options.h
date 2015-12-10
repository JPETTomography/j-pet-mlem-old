/// \page cmd_3d_tool 3D tools
/// \brief 3D helper utilities
///
/// Available utilities
/// -------------------
/// - \subpage cmd_3d_tool_crop
/// - \subpage cmd_3d_tool_psf

#pragma once

#include "cmdline.h"

namespace PET3D {
namespace Tool {

/// Adds \ref cmd_3d_tool_crop specific command line options.
void add_crop_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_tool_psf specific command line options.
void add_psf_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
void calculate_psf_options(cmdline::parser& cl, int argc = 1);

}  // Tool
}  // PET3D

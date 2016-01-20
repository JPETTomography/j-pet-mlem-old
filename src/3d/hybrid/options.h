/// \page cmd_3d_hybrid 3D hybrid
/// \brief 3D hybrid Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_3d_hybrid_matrix
/// - \subpage cmd_3d_hybrid_phantom
/// - \subpage cmd_3d_hybrid_reconstruction
/// - \subpage cmd_3d_hybrid_sensitivity
///
/// \sa \ref cmd_2d_barrel, \ref cmd_2d_strip

#pragma once

#include "common/options.h"

namespace PET3D {
namespace Hybrid {

/// Adds DetectorSet specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_hybrid_matrix specific command line options.
void add_matrix_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_hybrid_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_hybrid_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_hybrid_sensitivity specific command line options.
void add_sensitivity_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
void calculate_scanner_options(cmdline::parser& parser, int argc = 1);

/// Calculates all empty values from existing other parameters.
void calculate_phantom_options(cmdline::parser& cl, int argc = 1);

/// Calculates all empty values from existing other parameters.
void calculate_resonstruction_options(cmdline::parser& cl, int argc = 1);

/// Provides initialization list for creating detector.
#define __PET3D_LONGITUDINAL(...) __VA_ARGS__  // just pass-through
#define PET3D_LONGITUDINAL_SCANNER_CL(cl, ftype)                               \
  __PET3D_LONGITUDINAL(                                                        \
      Common::Convert<F, double>::cast(cl.get<std::vector<double>>("radius")), \
      Common::Convert<F, double>::cast(                                        \
          cl.get<std::vector<double>>("rotation")),                            \
      cl.get<std::vector<int>>("n-detectors"),                                 \
      cl.get<double>("w-detector"),                                            \
      cl.get<double>("h-detector"),                                            \
      cl.get<double>("d-detector"),                                            \
      cl.get<double>("fov-radius"))

enum Cmd { CmdReconstruction = 0, CmdPhantom, CmdPSF };
void calculate_cmd_options(cmdline::parser& cl, int argc, Cmd cmd);

}  // Hybrid
}  // PET3D

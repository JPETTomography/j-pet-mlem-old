/// \page cmd_3d_hybrid 3D Hybrid
/// \brief 3D Hybrid Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_3d_hybrid_matrix
/// - \subpage cmd_3d_hybrid_phantom
/// - \subpage cmd_3d_hybrid_reconstruction
/// - \subpage cmd_3d_hybrid_sensitivity
///
/// \sa \ref cmd_2d_barrel, \ref cmd_2d_strip
///
/// Workflow
/// --------
///
/// 1. Using 2D system matrix (recommended workflow):
///
/// \f[
///   \left.
///   \begin{array}{lll}
///     \mathit{scanner~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_matrix}
///                            &\!\!\!\!\rightarrow \mathit{2D~system~matrix}
///  \\ \mathit{phantom~desc.} &\!\!\!\!\rightarrow \mathtt{3d\_hybrid\_phantom}
///                            &\!\!\!\!\rightarrow \mathit{response}
///   \end{array}
///   \right\} \rightarrow \mathtt{3d\_hybrid\_reconstruction}
///            \rightarrow \mathit{reconstruction~image}
/// \f]
///
///    \note See \ref cmd_2d_barrel_matrix, \ref cmd_3d_hybrid_phantom and
///    \ref cmd_3d_hybrid_reconstruction for usage examples.
///
/// 2. Using 3D system matrix slice (invalid and **NOT** recommended workflow):
///
/// \f[
///   \left.
///   \begin{array}{lll}
///     \mathit{scanner~desc.} &\!\!\!\!\rightarrow \mathtt{3d\_hybrid\_matrix}
///                            &\!\!\!\!\rightarrow \mathit{system~matrix}
///  \\ \mathit{phantom~desc.} &\!\!\!\!\rightarrow \mathtt{3d\_hybrid\_phantom}
///                            &\!\!\!\!\rightarrow \mathit{response}
///   \end{array}
///   \right\} \rightarrow \mathtt{3d\_hybrid\_reconstruction}
///            \rightarrow \mathit{reconstruction~image}
/// \f]
///
/// 3. Using geometry description (**EXPERIMENTAL!**):
///
/// \f[
///   \left.
///   \begin{array}{lll}
///    \mathit{scanner~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_geometry}
///                           &\!\!\!\!\rightarrow \mathit{geometry~desc.}
/// \\ \mathit{phantom~desc.} &\!\!\!\!\rightarrow \mathtt{3d\_hybrid\_phantom}
///                           &\!\!\!\!\rightarrow \mathit{response}
///   \end{array}
///   \right\} \rightarrow \mathtt{3d\_hybrid\_reconstruction}
///            \rightarrow \mathit{reconstruction~image}
/// \f]

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

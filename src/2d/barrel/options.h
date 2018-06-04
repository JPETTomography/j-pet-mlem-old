/// \page cmd_2d_barrel 2D Barrel
/// \brief 2D Barrel Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_2d_barrel_matrix
/// - \subpage cmd_2d_barrel_phantom
/// - \subpage cmd_2d_barrel_reconstruction
/// - \subpage cmd_2d_barrel_geometry
/// - \subpage cmd_2d_barrel_lm_reconstruction
///
/// Workflow
/// --------
///
/// 1. Using system matrix:
///
/// \f[
///   \left.
///   \begin{array}{lll}
///     \mathit{scanner~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_matrix}
///                            &\!\!\!\!\rightarrow \mathit{system~matrix}
///  \\ \mathit{phantom~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_phantom}
///                            &\!\!\!\!\rightarrow \mathit{mean}
///   \end{array}
///   \right\} \rightarrow \mathtt{2d\_barrel\_reconstruction}
///            \rightarrow \mathit{reconstruction~image}
/// \f]
///
///    \note See \ref cmd_2d_barrel_matrix, \ref cmd_2d_barrel_phantom and
///    \ref cmd_2d_barrel_reconstruction for usage examples.
///
/// 2. Using LM and geometry description (**EXPERIMENTAL!**):
///
/// \f[
///   \left.
///   \begin{array}{lll}
///    \mathit{scanner~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_geometry}
///                           &\!\!\!\!\rightarrow \mathit{geometry~desc.}
/// \\ \mathit{phantom~desc.} &\!\!\!\!\rightarrow \mathtt{2d\_barrel\_phantom}
///                           &\!\!\!\!\rightarrow \mathit{response}
///   \end{array}
///   \right\} \rightarrow \mathtt{2d\_barrel\_lm\_reconstruction}
///            \rightarrow \mathit{reconstruction~image}
/// \f]

#pragma once

#include "common/options.h"

namespace PET2D {
namespace Barrel {

/// Adds scanner specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds custom config command line option.
void add_config_option(cmdline::parser& parser);

/// Adds pixel specific command line options.
void add_pixel_options(cmdline::parser& parser, bool required = false);

/// Adds probability model specific command line options.
void add_model_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_matrix specific command line options.
void add_matrix_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_lm_reconstruction specific command line options.
void add_lm_reconstruction_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
void calculate_scanner_options(cmdline::parser& parser,
                               int argc = 0,
                               bool calculate_pixel = true,
                               bool calculate_model = true);

/// Calculates all empty values from existing other parameters. (2nd version)
void calculate_scanner_options(cmdline::parser& cl,
                               int argc,
                               std::stringstream& assumed,
                               bool calculate_pixel = true,
                               bool calculate_model = true);

/// Provides initialization list for creating detector.
#define __PET2D_BARREL(...) __VA_ARGS__  // just pass-through
#define PET2D_BARREL_SCANNER_CL(cl, ftype)                                     \
  __PET2D_BARREL(                                                              \
      Common::Convert<F, double>::cast(cl.get<std::vector<double>>("radius")), \
      Common::Convert<F, double>::cast(                                        \
          cl.get<std::vector<double>>("rotation")),                            \
      cl.get<std::vector<int>>("n-detectors"),                                 \
      cl.get<double>("w-detector"),                                            \
      cl.get<double>("h-detector"),                                            \
      cl.get<double>("d-detector"),                                            \
      cl.get<double>("fov-radius"))

}  // Barrel
}  // PET2D

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
/// 2. Using LM and geometry description:
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

#include "cmdline.h"

namespace PET2D {
namespace Barrel {

/// Adds scanner specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds custom config command line option.
void add_config_option(cmdline::parser& parser);

/// Adds pixel specific command line options.
void add_pixel_options(cmdline::parser& parser);

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
void calculate_scanner_options(cmdline::parser& parser, int argc = 0);

/// Provides initialization list for creating detector.
#define __PET2D_BARREL(...) __VA_ARGS__  // just pass-through
#define PET2D_BARREL_SCANNER_CL(cl, ftype)               \
  __PET2D_BARREL({ (ftype)cl.get<double>("radius"),      \
                   (ftype)cl.get<double>("radius2"),     \
                   (ftype)cl.get<double>("radius3"),     \
                   (ftype)cl.get<double>("radius4") },   \
                 { (ftype)cl.get<double>("rotation"),    \
                   (ftype)cl.get<double>("rotation2"),   \
                   (ftype)cl.get<double>("rotation3"),   \
                   (ftype)cl.get<double>("rotation4") }, \
                 { cl.get<int>("n-detectors"),           \
                   cl.get<int>("n-detectors2"),          \
                   cl.get<int>("n-detectors3"),          \
                   cl.get<int>("n-detectors4") },        \
                 cl.get<double>("w-detector"),           \
                 cl.get<double>("h-detector"),           \
                 cl.get<double>("d-detector"))

}  // Barrel
}  // PET2D

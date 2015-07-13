/// \page cmd_3d_hybrid 3D hybrid
/// \brief 3D hybrid Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_3d_hybrid_matrix
// (disabled) - \subpage cmd_3d_hybrid_phantom
// (disabled) - \subpage cmd_3d_hybrid_reconstruction
///
/// \sa \ref cmd_2d_barrel, \ref cmd_2d_strip

#pragma once

#include "cmdline.h"

namespace PET3D {
namespace Hybrid {

/// Adds DetectorSet specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_hybrid_matrix specific command line options.
void add_matrix_options(cmdline::parser& parser);

// (disabled) Adds \ref cmd_3d_hybrid_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

// (disabled) Adds \ref cmd_3d_hybrid_reconstruction specific command line
// options.
void add_reconstruction_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
void calculate_scanner_options(cmdline::parser& parser);

void set_small_barrel_options(cmdline::parser& parser);
void set_big_barrel_options(cmdline::parser& parser);

/// Provides initialization list for creating detector.
#define __PET3D_LONGITUDINAL(...) __VA_ARGS__  // just pass-through
#define PET3D_LONGITUDINAL_SCANNER_CL(cl, ftype)               \
  __PET3D_LONGITUDINAL({ (ftype)cl.get<double>("radius"),      \
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

}  // Hybrid
}  // PET3D

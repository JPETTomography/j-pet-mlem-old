/// \page cmd_3d_longitudinal 3D longitudinal
/// \brief 3D longitudinal Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_3d_longitudinal_matrix
/// - \subpage cmd_3d_longitudinal_phantom
/// - \subpage cmd_3d_longitudinal_reconstruction
///
/// \sa \ref cmd_2d_barrel, \ref cmd_2d_strip

#pragma once

#include "cmdline.h"

namespace PET3D {
namespace Longitudinal {

/// Adds DetectorSet specific command line options.
void add_detector_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_longitudinal_matrix specific command line options.
void add_matrix_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_longitudinal_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// Adds \ref cmd_3d_longitudinal_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
void calculate_detector_options(cmdline::parser& parser);

void set_small_barrel_options(cmdline::parser& parser);
void set_big_barrel_options(cmdline::parser& parser);

/// Provides initialization list for creating detector.
#define __PET3D_LONGITUDINAL(...) __VA_ARGS__  // just pass-through
#define PET3D_LONGITUDINAL_DETECTOR_CL(cl, ftype)              \
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

}  // Strip
}  // PET2D

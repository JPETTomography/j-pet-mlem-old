/// \page cmd_2d_barrel 2D Barrel
/// \brief 2D Barrel Monte-Carlo & reconstruction tools
///
/// Available tools
/// ---------------
/// - \subpage cmd_2d_barrel_matrix
/// - \subpage cmd_2d_barrel_phantom
/// - \subpage cmd_2d_barrel_reconstruction

#pragma once

#include "cmdline.h"

namespace PET2D {
namespace Barrel {

/// Adds scanner specific command line options.
void add_scanner_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_matrix specific command line options.
void add_matrix_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_phantom specific command line options.
void add_phantom_options(cmdline::parser& parser);

/// Adds \ref cmd_2d_barrel_reconstruction specific command line options.
void add_reconstruction_options(cmdline::parser& parser);

/// Calculates all empty values from existing other parameters.
//void calculate_scanner_options(cmdline::parser& parser);

void set_small_barrel_options(cmdline::parser& parser);
void set_big_barrel_options(cmdline::parser& parser);

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

}  // Strip
}  // PET2D

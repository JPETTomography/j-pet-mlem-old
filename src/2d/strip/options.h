#pragma once

#include "cmdline.h"

namespace PET2D {
namespace Strip {

/// adds reconstruction specific command line options
void add_reconstruction_options(cmdline::parser& parser);

/// adds phantom specific command line options
void add_phantom_options(cmdline::parser& parser);

/// calculates all empty values from existing other parameters
void calculate_detector_options(cmdline::parser& parser);

/// provides initialization list for creating detector
#define CL_DETECTOR_PARAMETERS(cl)                          \
  cl.get<double>("r-distance"), cl.get<double>("s-length"), \
      cl.get<int>("n-y-pixels"), cl.get<int>("n-z-pixels"), \
      cl.get<double>("p-size"), cl.get<double>("p-size"),   \
      cl.get<double>("s-z"), cl.get<double>("s-dl")

}  // Strip
}  // PET2D

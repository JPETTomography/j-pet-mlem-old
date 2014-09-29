#ifndef OPTIONS_H
#define OPTIONS_H

#include "cmdline.h"

void set_detector_options(cmdline::parser& parser);
void set_options_for_reconstruction(cmdline::parser& parser);
void set_options_for_phantom(cmdline::parser& parser);

#endif  // OPTIONS_H

#include "options.h"

#include "2d/barrel/options.h"

namespace PET3D {
namespace Hybrid {

void add_detector_options(cmdline::parser& cl) {
  PET2D::Barrel::add_detector_options(cl);
}

void add_matrix_options(cmdline::parser& cl) {
  PET2D::Barrel::add_matrix_options(cl);
}

void add_phantom_options(cmdline::parser& cl) {
  PET2D::Barrel::add_phantom_options(cl);
}

void add_reconstruction_options(cmdline::parser& cl) {
  PET2D::Barrel::add_reconstruction_options(cl);
}

void calculate_detector_options(cmdline::parser& cl) {
  PET2D::Barrel::calculate_detector_options(cl);
}
}
}

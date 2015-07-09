#include "options.h"

#include "2d/barrel/options.h"

namespace PET3D {
namespace Hybrid {

void add_scanner_options(cmdline::parser& cl) {
  PET2D::Barrel::add_scanner_options(cl);
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

void calculate_scanner_options(cmdline::parser& cl) {
  PET2D::Barrel::calculate_scanner_options(cl);
}

void set_small_barrel_options(cmdline::parser& parser) {
  PET2D::Barrel::set_small_barrel_options(parser);
  parser.get<double>("length") = 0.3;
}

void set_big_barrel_options(cmdline::parser& parser) {
  PET2D::Barrel::set_big_barrel_options(parser);
  parser.get<double>("length") = 0.5;
}
}
}

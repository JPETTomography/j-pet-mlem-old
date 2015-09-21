#include "options.h"

#include "2d/barrel/options.h"
#include "util/cmdline_types.h"

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
  cl.add<std::string>("geometry", 0, "geometry information", true);
  cl.add<double>("s-z", 0, "TOF sigma along z axis", false, 0.015);
  cl.add<double>("length", 0, "length of the detector", false, 0.3);

  PET2D::Barrel::add_matrix_options(cl);
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);

  cl.footer("response ...");
}

void add_sensitivity_options(cmdline::parser& cl) {
  PET2D::Barrel::add_scanner_options(cl);
  cl.add<int>("z-plane", 0, "z plane trianguler cut", false);
  cl.add<int>("y-plane", 0, "y plane cut", false);
  cl.add<int>(
      "n-pixels", 'n', "number of pixels in x and y  directions", false, 80);
  cl.add<int>("n-planes", 0, "number pf z planes", false, 80);
  cl.add<double>("s-pixel", 'p', "voxel size", false, 0.005);
  cl.add<int>("n-emissions", 'e', "number of emission", false, 0);
  cl.add<cmdline::path>(
      "output", 'o', "output files template", false, "out.bin");
  cl.add<int>("n-threads", 'T', "number of threads", false);
}

void calculate_scanner_options(cmdline::parser& cl, int argc) {
  PET2D::Barrel::calculate_scanner_options(cl, argc);
}

}  // Hybrid
}  // PET3D

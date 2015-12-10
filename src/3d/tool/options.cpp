#include "options.h"

#include "common/options.h"
#include "3d/hybrid/options.h"
#include "2d/barrel/options.h"

#include "util/cmdline_types.h"

namespace PET3D {
namespace Tool {

void add_psf_options(cmdline::parser& cl) {
  PET2D::Barrel::add_pixel_options(cl, true);
  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  Common::add_openmp_options(cl);
  cl.footer("image ...");
}

void add_crop_options(cmdline::parser& cl) {
  add_psf_options(cl);
  cl.add<int>("crop", 'c', "numer of pixel to crop the output image", true);
  cl.add<int>("crop-x", 'x', "crop origin pixel x", false);
  cl.add<int>("crop-y", 'y', "crop origin pixel y", false);
  cl.add<int>("crop-z", 'z', "crop origin pixel z", false);
  cl.add<cmdline::path>("output", 'o', "output image", true);
}

void calculate_psf_options(cmdline::parser& cl, int argc) {
  PET3D::Hybrid::calculate_cmd_options(cl, argc, Hybrid::CmdPSF);
}

}  // Tool
}  // PET3D

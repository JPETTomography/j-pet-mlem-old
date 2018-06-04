#include "options.h"

#include "common/options.h"
#include "3d/hybrid/options.h"
#include "2d/barrel/options.h"

#include "util/cmdline_types.h"

#include "util/random.h"  // util::random::tausworthe::seed_type

namespace PET3D {
namespace Tool {

void add_psf_options(cmdline::parser& cl) {
  PET2D::Barrel::add_pixel_options(cl, true);
  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  cl.add<int>("padding", 0, "pad given pixels seeking maximum", false, 0);
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

void add_convert_options(cmdline::parser& cl) {
  cl.footer("means");

  PET2D::Barrel::add_config_option(cl);
  PET2D::Barrel::add_scanner_options(cl);

  cl.add<double>("s-dl", 0, "TOF sigma delta-l", cmdline::alwayssave, 0.06);
  cl.add<double>(
      "s-z", 0, "TOF sigma along z axis", cmdline::alwayssave, 0.015);

  cl.add<double>("length", 'l', "scintillator length", false, 1.0);

  // printing & stats params
  cl.add("verbose", 'v', "prints the iterations information on std::out");
  cl.add<util::random::tausworthe::seed_type>(
      "seed", 'S', "random number generator seed", cmdline::dontsave);

  // conversion types
  cl.add("warsaw", 0, "Warsaw data format", cmdline::dontsave);
  cl.add("tbednarski", 0, "T.Bednarski data format", cmdline::dontsave);

  // Warsaw data specific
  cl.add<int>("only",
              0,
              "limit event (scattering) type, eg. 1 for non-scattered",
              false,
              0);

  // T.Bednarski data specific: for debugging precision
  cl.add("relative-time",
         0,
         "output relative time instead z_u, z_d, dl",
         cmdline::dontsave);
}

void calculate_psf_options(cmdline::parser& cl, int argc) {
  PET3D::Hybrid::calculate_cmd_options(cl, argc, Hybrid::CmdPSF);
}

}  // Tool
}  // PET3D

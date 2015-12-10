#include "options.h"

#include "2d/barrel/options.h"
#include "2d/strip/options.h"
#include "common/options.h"
#include "util/cmdline_types.h"

namespace PET3D {
namespace Hybrid {

void add_scanner_options(cmdline::parser& cl) {
  PET2D::Barrel::add_scanner_options(cl);
}

void add_matrix_options(cmdline::parser& cl) {
  PET2D::Barrel::add_matrix_options(cl);
  cl.add<double>("z-position", 'z', "position of the z plane", false, 0);
  cl.add<double>("length", 'L', "length of the detector", false, 2);
}

void add_phantom_options(cmdline::parser& cl) {
  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  cl.add<double>("length", 'L', "length of the detector", false, 2);
  cl.add<double>("s-z", 0, "TOF sigma along z axis", false, 0.015);
  cl.add<double>("z-position", 'z', "position of the z plane", false, 0);
  PET2D::Barrel::add_phantom_options(cl);
}

void add_reconstruction_options(cmdline::parser& cl) {
  cl.add<std::string>("geometry", 0, "geometry information", true);
  cl.add<cmdline::path>("system", 0, "system matrix file", false);
  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  cl.add<double>("length", 'L', "length of the detector", false, 2);
  cl.add<double>("s-z", 0, "TOF sigma along z axis", false, 0.015);
  cl.add("sens-to-one", 0, "sets sensitivity to one", false);

  cl.add<int>(
      "crop-pixels", 0, "numer of pixel to crop the output image", false);
  cl.add<int>("crop-x", 0, "crop origin pixel x", false);
  cl.add<int>("crop-y", 0, "crop origin pixel y", false);
  cl.add<int>("crop-z", 0, "crop origin pixel z", false);

  PET2D::Barrel::add_matrix_options(cl);
  cl.add<int>("blocks", 'i', "number of iteration blocks", false, 0);
  cl.add<int>("iterations", 'I', "number of iterations (per block)", false, 1);
  cl.add<double>("z-position", 'z', "position of the z plane", false, 0);
  cl.add<cmdline::path>("rho", 0, "start rho (eg. existing iteration)", false);
  cl.footer("response ...");
}

void add_sensitivity_options(cmdline::parser& cl) {

  PET2D::Barrel::add_matrix_options(cl);

  cl.add<int>("z-plane", 0, "z plane trianguler cut", false);
  cl.add<int>("y-plane", 0, "y plane cut", false);
  cl.add<int>("n-planes", 0, "number of z planes", false, 80);
  cl.add<double>("length", 'L', "length of the detector", false, 2);
  cl.add("int", 'i', "values are integers");

  cl.footer("image ...");
}

void add_psf_options(cmdline::parser& cl) {

  PET2D::Barrel::add_pixel_options(cl, true);
  cl.add<int>("n-planes", 0, "number of voxels in z direction", false, 0);
  cl.add<double>("z-left", 0, "left extent in z direction", false, 0);
  Common::add_openmp_options(cl);
}

void calculate_scanner_options(cmdline::parser& cl, int argc) {
  PET2D::Barrel::calculate_scanner_options(cl, argc);
}

enum Cmd { CmdReconstruction = 0, CmdPhantom, CmdPSF };

static void calculate_cmd_options(cmdline::parser& cl, int argc, Cmd cmd) {
  std::stringstream assumed;
  auto calculate_pixel = cmd != CmdPhantom || cl.exist("n-pixels");
  if (cmd != CmdPSF) {
    PET2D::Barrel::calculate_scanner_options(
        cl, argc, assumed, calculate_pixel);
  }

  auto& n_pixels = cl.get<int>("n-pixels");
  auto& n_planes = cl.get<int>("n-planes");
  auto& s_pixel = cl.get<double>("s-pixel");
  auto& z_left = cl.get<double>("z-left");

  if (calculate_pixel) {
    if (!cl.exist("n-planes")) {
      n_planes = n_pixels;
      assumed << "--n-planes=" << n_planes << std::endl;
    }
    if (!cl.exist("z-left")) {
      z_left = -n_planes * s_pixel / 2;
      assumed << "--z-left=" << z_left << std::endl;
    }
  }

  if ((cmd == CmdPSF || cl.exist("verbose")) && assumed.str().size()) {
    std::cerr << "assumed:" << std::endl << assumed.str();
  }
}

void calculate_phantom_options(cmdline::parser& cl, int argc) {
  calculate_cmd_options(cl, argc, CmdPhantom);
}

void calculate_resonstruction_options(cmdline::parser& cl, int argc) {
  calculate_cmd_options(cl, argc, CmdReconstruction);
}

void calculate_psf_options(cmdline::parser& cl, int argc) {
  calculate_cmd_options(cl, argc, CmdPSF);
}

}  // Hybrid
}  // PET3D

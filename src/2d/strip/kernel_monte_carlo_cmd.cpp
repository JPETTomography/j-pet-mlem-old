#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "common/types.h"

#include "2d/strip/kernel_monte_carlo.h"
#include "2d/geometry/pixel_map.h"
#include "2d/geometry/pixel.h"

#include "2d/strip/options.h"
#include "2d/strip/response.h"
#include "2d/strip/reconstruction.h"
#include "2d/strip/gaussian_kernel.h"

using Pixel = PET2D::Pixel<int>;
using Output = PET2D::PixelMap<Pixel, F>;

using Kernel = PET2D::Strip::GaussianKernel<F>;
using Reconstruction = PET2D::Strip::Reconstruction<F, Kernel>;
using Scanner = PET2D::Strip::Scanner<F, S>;

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  PET2D::Strip::add_reconstruction_options(cl);
  cl.add<cmdline::path>("rho", 0, "start rho (eg. existing iteration)", false);
  cl.parse_check(argc, argv);
  PET2D::Strip::calculate_scanner_options(cl, argc);
  if (argc == 1) {
    strip_integral();
    strip_integral_theta();
    strip_integral_event();
    strip_integral_theta_event();
    strip_gauss_kernel();
    strip_gauss_kernel_integral();
  } else {

    Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));
    Reconstruction reconstruction(scanner);

    auto verbose = cl.count("verbose");
    if (verbose) {
      std::cout << "# image: " << scanner.n_y_pixels << "x"
                << scanner.n_z_pixels << std::endl;
    }

    auto is_3d = cl.exist("3d-response");

    for (auto& fn : cl.rest()) {
      if (cmdline::path(fn).ext() == ".txt") {
#if USE_FAST_TEXT_PARSER
        reconstruction.fast_load_txt_events(fn.c_str(), is_3d);
#else
        if (is_3d) {
          throw("3D input not supported in this build");
        }
        std::ifstream in_responses(fn);
        ENSURE_IS_OPEN(in_responses, "phantom responses", fn);
        reconstruction << in_responses;
#endif
      } else {
        if (is_3d) {
          throw("3D input must have .txt extension");
        }
        util::ibstream in_responses(fn);
        ENSURE_IS_OPEN(in_responses, "phantom responses", fn);
        reconstruction << in_responses;
      }

      if (verbose) {
        std::cerr << "# read " << reconstruction.responses.size()
                  << " responsess from " << fn << std::endl;
      }
    }

    if (cl.exist("rho")) {
      auto rho_name = cl.get<cmdline::path>("rho");
      auto rho_base_name = rho_name.wo_ext();
      auto rho_ext = rho_name.ext();
      auto rho_txt = rho_ext == ".txt";
      if (rho_txt) {
        std::ifstream txt(rho_name);
        txt >> reconstruction.rho;
      } else {
        util::ibstream bin(rho_name);
        bin >> reconstruction.rho;
      }
    }

    auto output_name = cl.get<cmdline::path>("output");
    auto output_base_name = output_name.wo_ext();
    auto output_ext = output_name.ext();
    auto output_txt = output_ext == ".txt";

    // FIXME: implement me!
    (void)output_txt;  // stop unused warning
  }
}

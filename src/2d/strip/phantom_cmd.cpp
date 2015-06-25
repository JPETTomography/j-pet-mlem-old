/// \page cmd_2d_strip_phantom 2d_strip_phantom
/// \brief 2D Strip PET phantom tool
///
/// Simulates scanner response for given virtual phantom and produces mean file
/// for \ref cmd_2d_strip_reconstruction.
///
/// Example phantom descriptions
/// ----------------------------
/// - Shepp like phantom
///
///   \verbinclude phantom/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantom/s_shepp_small
///
/// Authors
/// -------
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/strip/phantom_cmd.txt
///
/// \sa \ref cmd_2d_strip_reconstruction

#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <random>

#if SSE_FLUSH
#include <xmmintrin.h>
#endif

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "options.h"

#include "2d/strip/phantom.h"
#include "2d/strip/scanner.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

using F = float;
using S = short;
using Hit = int;

using Scanner = PET2D::Strip::Scanner<F, S>;
using Phantom = PET2D::Strip::Phantom<F, S>;
using Ellipse = PET2D::Ellipse<F>;
using PhantomRegion = PET2D::Strip::PhantomRegion<F>;

const double RADIAN = M_PI / 180;

int main(int argc, char* argv[]) {

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  try {
    cmdline::parser cl;
    PET2D::Strip::add_phantom_options(cl);
    cl.parse_check(argc, argv);
    PET2D::Strip::calculate_scanner_options(cl);

    if (!cl.rest().size()) {
      throw(
          "at least one input phantom description file expected, "
          "consult --help");
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    auto emissions = cl.get<int>("emissions");
    auto verbose = cl.exist("verbose");

    std::vector<PhantomRegion> ellipse_list;

    Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));

    if (verbose) {
      std::cerr << "size: " << scanner.n_z_pixels << "x" << scanner.n_y_pixels
                << std::endl;
    }

    for (auto& fn : cl.rest()) {
      std::ifstream infile(fn);
      std::string line;
      while (std::getline(infile, line)) {
        std::istringstream iss(line);
        std::string type;
        iss >> type;
        if (type == "ellipse") {
          double x, y, a, b, angle, acceptance;

          // on error
          if (!(iss >> x >> y >> a >> b >> angle >> acceptance))
            break;

          Ellipse el(x, y, a, b, angle * RADIAN);

          if (verbose) {
            std::cout << "ellipse: " << el.center.x << " " << el.center.y << " "
                      << el.a << " " << el.b << " " << el.angle << " " << el.A
                      << " " << el.B << " " << el.C << std::endl;
          }

          PhantomRegion region(el, acceptance);
          ellipse_list.push_back(region);
        } else {
          std::cerr << "unknow phantom type" << std::endl;
          exit(-1);
        }
      }
    }

    Phantom phantom(ellipse_list);

    if (verbose) {
      std::cerr << "scanner: " << scanner.size_y << " " << scanner.tl_y_half_h
                << std::endl;
    }

    Common::PhantomMonteCarlo<Phantom, Scanner> monte_carlo(phantom, scanner);

    typename Phantom::RNG rng;
    Common::AlwaysAccept<F> model;

    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    std::ofstream out_wo_error(output_base_name + "_geom_only" + ext);
    monte_carlo.out_wo_error = out_wo_error;

    std::ofstream out_w_error(output);
    monte_carlo.out_w_error = out_w_error;

    std::ofstream out_exact_events(output_base_name + "_exact_events" + ext);
    monte_carlo.out_exact_events = out_exact_events;

    std::ofstream out_full_response(output_base_name + "_full_response" + ext);
    monte_carlo.out_full_response = out_full_response;

    monte_carlo.generate(rng, model, emissions);
    monte_carlo.write_out(rng);

    if (verbose) {
      std::cerr << "detected: " << monte_carlo.n_events_detected() << " events"
                << std::endl;
    }

    std::ofstream cfg(output_base_name + ".cfg");
    cfg << cl;

  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }

  return 0;
}

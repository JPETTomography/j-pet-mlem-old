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
#include "util/random.h"
#include "util/backtrace.h"
#include "options.h"

#include "2d/geometry/phantom.h"
#include "2d/strip/scanner.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

using F = float;
using S = short;
using Hit = int;

using RNG = util::random::tausworthe;
using Scanner = PET2D::Strip::Scanner<F, S>;
using Phantom = PET2D::Phantom<RNG, F>;
using Ellipse = PET2D::Ellipse<F>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;

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

    auto n_emissions = cl.get<int>("n-emissions");
    auto verbose = cl.exist("verbose");

    Scanner scanner(PET2D_STRIP_SCANNER_CL(cl));

    if (verbose) {
      std::cerr << "size: " << scanner.n_z_pixels << "x" << scanner.n_y_pixels
                << std::endl;
    }

    Phantom phantom;
    for (auto& fn : cl.rest()) {
      std::ifstream infile(fn);
      phantom.read_from_stream(infile);
    }

    phantom.calculate_cdf();

    if (verbose) {
      std::cerr << "scanner: " << scanner.size_y << " " << scanner.tl_y_half_h
                << std::endl;
    }

    MonteCarlo monte_carlo(phantom, scanner);

    typename Phantom::RNG rng;
    Common::AlwaysAccept<F> model;

    if (cl.exist("output")) {
      auto output = cl.get<cmdline::path>("output");
      auto output_base_name = output.wo_ext();
      auto ext = output.ext();

      std::ofstream out_wo_error(output_base_name + "_geom_only" + ext);
      std::ofstream out_w_error(output);
      std::ofstream out_exact_events(output_base_name + "_exact_events" + ext);
      std::ofstream out_full_response(output_base_name + "_full_response" +
                                      ext);

      monte_carlo(
          rng,
          model,
          n_emissions,
          [](Phantom::Event&) {},
          [&](const typename MonteCarlo::Event& event,
              const typename MonteCarlo::FullResponse& full_response) {
            out_exact_events << event << "\n";
            out_full_response << full_response << "\n";
            out_wo_error << scanner.response_wo_error(full_response) << "\n";
            out_w_error << scanner.response_w_error(rng, full_response) << "\n";
          });
      if (verbose) {
        std::cerr << "detected: " << monte_carlo.n_events_detected()
                  << " events" << std::endl;
      }
      std::ofstream cfg(output_base_name + ".cfg");
      cfg << cl;
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
    util::print_backtrace(std::cerr);
  }

  return 0;
}

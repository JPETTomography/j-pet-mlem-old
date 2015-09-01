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
///   \verbinclude phantoms/s_shepp
///
/// - Small Shepp like phantom
///
///   \verbinclude phantoms/s_shepp_scaled
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
#include "util/cmdline_hooks.h"
#include "util/random.h"
#include "util/backtrace.h"
#include "util/progress.h"
#include "util/png_writer.h"
#include "options.h"

#include "2d/geometry/phantom.h"
#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/pixel_map.h"
#include "2d/strip/scanner.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"
#include "common/types.h"

#if _OPENMP
#include <omp.h>
#endif

using RNG = util::random::tausworthe;
using Scanner = PET2D::Strip::Scanner<F, S>;
using Phantom = PET2D::Phantom<RNG, F>;
using Ellipse = PET2D::Ellipse<F>;
using Image = PET2D::PixelMap<PET2D::Pixel<S>, int>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner, Image>;
using Event = MonteCarlo::Event;
using FullResponse = MonteCarlo::FullResponse;

int main(int argc, char* argv[]) {
  CMDLINE_TRY

#if SSE_FLUSH
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
#endif

  cmdline::parser cl;
  PET2D::Strip::add_phantom_options(cl);
  cl.parse_check(argc, argv);
  PET2D::Strip::calculate_scanner_options(cl, argc);

  if (!cl.rest().size()) {
    if (argc == 1) {
      std::cerr << cl.usage();
      exit(0);
    } else {
      throw("at least one input phantom description expected, consult --help");
    }
  }

#if _OPENMP
  if (cl.exist("n-threads")) {
    omp_set_num_threads(cl.get<int>("n-threads"));
  }
#endif

  auto n_emissions = cl.get<int>("n-emissions");
  auto verbose = cl.count("verbose");

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

  RNG rng;
  Common::AlwaysAccept<F> model;

  if (!cl.exist("output")) {
    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    std::ofstream out_wo_error(output_base_name + "_wo_error" + ext);
    std::ofstream out_w_error(output);
    std::ofstream out_exact_events(output_base_name + "_events" + ext);
    std::ofstream out_full_response(output_base_name + "_full_response" + ext);

    auto n_z_pixels = cl.get<int>("n-z-pixels");
    auto n_y_pixels = cl.get<int>("n-y-pixels");
    auto s_pixel = cl.get<double>("s-pixel");
    PET2D::PixelGrid<F, S> pixel_grid(
        n_z_pixels,
        n_y_pixels,
        s_pixel,
        PET2D::Point<F>(-s_pixel * n_z_pixels / 2, -s_pixel * n_y_pixels / 2));

    Image image_emitted(n_z_pixels, n_y_pixels);
    Image image_detected_exact(n_z_pixels, n_y_pixels);
    Image image_detected_w_error(n_z_pixels, n_y_pixels);

    util::progress progress(verbose, n_emissions, 10000);
    monte_carlo(
        rng,
        model,
        n_emissions,
        [&](const Event& event) {
          auto pixel = pixel_grid.pixel_at(event.center);
          if (pixel_grid.contains(pixel)) {
            image_emitted[pixel]++;
          }
        },
        [&](const Event& event, const FullResponse& full_response) {
          out_exact_events << event << "\n";
          out_full_response << full_response << "\n";
          out_wo_error << scanner.response_wo_error(full_response) << "\n";
          auto response_w_error = scanner.response_w_error(rng, full_response);
          out_w_error << response_w_error << "\n";
          {
            auto pixel = pixel_grid.pixel_at(event.center);
            if (pixel_grid.contains(pixel)) {
              image_detected_exact[pixel]++;
            }
          }
          {
            auto event = scanner.from_projection_space_tan(response_w_error);
            auto pixel = pixel_grid.pixel_at(PET2D::Point<F>(event.z, event.y));
            if (pixel_grid.contains(pixel)) {
              image_detected_w_error[pixel]++;
            }
          }
        },
        progress);
    if (verbose) {
      std::cerr << " emitted: " << monte_carlo.n_events_emitted() << " events"
                << std::endl
                << "detected: " << monte_carlo.n_events_detected() << " events"
                << std::endl;
    }

    std::ofstream cfg(output_base_name + ".cfg");
    cfg << cl;

    util::png_writer png_wo_error(output_base_name + "_wo_error.png");
    png_wo_error << image_detected_exact;
    util::png_writer png_emitted(output_base_name + "_emitted.png");
    png_emitted << image_emitted;
    util::png_writer png_w_error(output_base_name + ".png");
    png_w_error << image_detected_w_error;
  }

  CMDLINE_CATCH
}


#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/backtrace.h"
#include "util/progress.h"

#include "2d/geometry/pixel_grid.h"
#include "2d/geometry/pixel_map.h"
#include "2d/geometry/phantom.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"
#include "2d/toy/gauss_scanner.h"
#include "common/types.h"

using RNG = util::random::tausworthe;
using Scanner = PET2D::Toy::GaussScanner<F>;
using Phantom = PET2D::Phantom<RNG, F>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;
using Event = MonteCarlo::Event;
using Image = PET2D::PixelMap<PET2D::Pixel<S>, Hit>;
using FullResponse = MonteCarlo::FullResponse;

int main(int argc, char* argv[]) {

  CMDLINE_TRY
  cmdline::parser cl;

  cl.add<double>("length-x", 'l', "length in x", false, 1);
  cl.add<double>("length-y", '\0', "length in y", false, 1);

  cl.add<double>("s-pixel", 'p', "pixel size", false, 0.001);
  cl.add<int>("n-pixels", 'n', "number of pixels", cmdline::dontsave, 0);
  cl.add<int>("n-z-pixels", 0, "number of z pixels", false);
  cl.add<int>("n-y-pixels", 0, "number of y pixels", false);
  cl.add<double>(
      "s-x", 0, "TOF sigma along x axis", cmdline::alwayssave, 0.012);
  cl.add<double>("s-y", 0, "TOF sigma along y axis", cmdline::alwayssave, 0.04);
  cl.add<int>("n-emissions", 'e', "number of emissions", false, 0);
  cl.add<double>("scale", '\0', "scale factor", false, 1);
  cl.add("verbose", 'v', "print progress information (-v)");
  cl.add("detected", '\0', "count only detected events");
  cl.add<cmdline::path>("output",
                        'o',
                        "output files prefix (png)",
                        false,
                        cmdline::path(),
                        cmdline::not_from_file);
  cl.add("no-responses", 0, "Do not emit responses", cmdline::dontsave);
#if _OPENMP
  cl.add<int>("n-threads", 'T', "number of OpenMP threads", cmdline::dontsave);
#else
  (void)cl;  // mark cl as unsued when not using OpenMP
#endif

  cl.parse_check(argc, argv);

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

  Phantom phantom(cl.get<double>("scale"));
  for (auto& fn : cl.rest()) {
    std::ifstream in_phantom(fn);
    phantom << in_phantom;
  }

  phantom.calculate_cdf();

  auto n_emissions = cl.get<int>("n-emissions");
  auto verbose = cl.count("verbose");

  Scanner scanner(cl.get<double>("s-x"), cl.get<double>("s-y"));
  if (verbose)
    std::cout << "sigma : " << cl.get<double>("s-x") << " "
              << cl.get<double>("s-y") << "\n";

  MonteCarlo monte_carlo(phantom, scanner);

  RNG rng;
  Common::AlwaysAccept<F> model;

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();
  bool no_responses = cl.exist("no-responses");

  std::ofstream out_wo_error, out_w_error, out_exact_events, out_full_response;
  if (output_base_name.length() && !no_responses) {
    out_wo_error.open(output_base_name + "_wo_error" + ext);
    out_w_error.open(output);
    out_exact_events.open(output_base_name + "_events" + ext);
    out_full_response.open(output_base_name + "_full_response" + ext);
  } else {
    no_responses = true;
  }

  auto n_z_pixels = cl.get<int>("n-pixels");
  auto n_y_pixels = cl.get<int>("n-pixels");
  auto s_pixel = cl.get<double>("s-pixel");

  PET2D::PixelGrid<F, S> pixel_grid(
      n_z_pixels,
      n_y_pixels,
      s_pixel,
      PET2D::Point<F>(-s_pixel * n_z_pixels / 2, -s_pixel * n_y_pixels / 2));

  Image image_emitted(n_z_pixels, n_y_pixels);
  Image image_detected_exact(n_z_pixels, n_y_pixels);
  Image image_detected_w_error(n_z_pixels, n_y_pixels);

  bool only_detected = cl.exist("detected");

  util::progress progress(verbose, n_emissions, 10000);
  monte_carlo(
      rng,
      model,
      n_emissions,
      [&](const Event& event) {
        auto pixel = pixel_grid.pixel_at(event.origin);
        if (pixel_grid.contains(pixel)) {
          image_emitted[pixel]++;
        }
      },
      [&](const Event& event, const FullResponse& full_response) {
        auto response_w_error = scanner.response_w_error(rng, full_response);
        if (!no_responses) {
          out_exact_events << event << "\n";
          out_full_response << full_response << "\n";
          out_wo_error << scanner.response_wo_error(full_response) << "\n";
          out_w_error << response_w_error << "\n";
        }
        {
          auto pixel = pixel_grid.pixel_at(event.origin);
          if (pixel_grid.contains(pixel)) {
            image_detected_exact[pixel]++;
          }
        }
        {

          auto pixel = pixel_grid.pixel_at(
              PET2D::Point<F>(response_w_error.x, response_w_error.y));
          if (pixel_grid.contains(pixel)) {
            image_detected_w_error[pixel]++;
          }
        }
      },
      progress,
      only_detected);
  if (verbose) {
    std::cerr << " emitted: " << monte_carlo.n_events_emitted() << " events"
              << std::endl
              << "detected: " << monte_carlo.n_events_detected() << " events"
              << std::endl;
  }

  if (output_base_name.length()) {
    std::ofstream cfg(output_base_name + ".cfg");
    cfg << cl;

    // RAW + NRRD
    util::obstream bin_wo_error(output_base_name + "_wo_error");
    util::nrrd_writer nrrd_wo_error(output_base_name + "_wo_error.nrrd",
                                    output_base_name + "_wo_error");
    bin_wo_error << image_detected_exact;
    nrrd_wo_error << image_detected_exact;
    util::obstream bin_emitted(output_base_name + "_emitted");
    util::nrrd_writer nrrd_emitted(output_base_name + "_emitted.nrrd",
                                   output_base_name + "_emitted");
    bin_emitted << image_emitted;
    nrrd_emitted << image_emitted;
    util::obstream bin_w_error(output_base_name + "_w_error");
    util::nrrd_writer nrrd_w_error(output_base_name + "_w_error.nrrd",
                                   output_base_name + "_w_error");
    bin_w_error << image_detected_w_error;
    nrrd_w_error << image_detected_w_error;
  }

  // PNG
  util::png_writer png_wo_error(output_base_name + "_wo_error.png");
  png_wo_error << image_detected_exact;
  util::png_writer png_emitted(output_base_name + "_emitted.png");
  png_emitted << image_emitted;
  util::png_writer png_w_error(output_base_name + "_w_error.png");
  png_w_error << image_detected_w_error;

  CMDLINE_CATCH
}

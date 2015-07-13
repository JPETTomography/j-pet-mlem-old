/// \page cmd_2d_barrel_phantom 2d_barrel_phantom
/// \brief 2D Barrel PET phantom generation tool
///
/// Simulates detector response for given virtual phantom and produces mean file
/// for \ref cmd_2d_barrel_reconstruction.
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
/// - Piotr Bialas    <piotr.bialas@uj.edu.pl>
/// - Jakub Kowal     <jakub.kowal@uj.edu.pl>
/// - Adam Strzelecki <adam.strzelecki@uj.edu.pl>
///
/// Usage
/// -----
/// \verbinclude src/2d/barrel/phantom_cmd.txt
///
/// \sa \ref cmd_2d_barrel_matrix, \ref cmd_2d_barrel_reconstruction

#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"

#include "2d/geometry/point.h"
#include "2d/barrel/scanner_builder.h"
#include "ring_scanner.h"
#include "generic_scanner.h"
#include "circle_detector.h"
#include "triangle_detector.h"
#include "polygonal_detector.h"

#include "util/png_writer.h"
#include "util/progress.h"
#include "util/json.h"
#include "util/random.h"
#include "options.h"

#include "2d/geometry/phantom.h"
#include "common/model.h"
#include "common/phantom_monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

using F = float;
using S = short;
using Hit = int;

using RNG = util::random::tausworthe;

using Pixel = PET2D::Pixel<S>;
using Point = PET2D::Point<F>;
using Event = PET2D::Event<F>;

template <class DetectorClass>
using Scanner = PET2D::Barrel::GenericScanner<DetectorClass, S>;
template <class DetectorClass>
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<DetectorClass>;

// all available detector shapes
using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;
using CircleScanner = Scanner<PET2D::Barrel::CircleDetector<F>>;
using TriangleScanner = Scanner<PET2D::Barrel::TriangleDetector<F>>;
using HexagonalScanner = Scanner<PET2D::Barrel::PolygonalDetector<6, F>>;

using Ellipse = PET2D::Ellipse<F>;
using Phantom = PET2D::Phantom<RNG, F>;

template <class DetectorClass, class PhantomClass, class ModelClass>
void run(cmdline::parser& cl, PhantomClass& phantom, ModelClass& model);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    PET2D::Barrel::add_phantom_options(cl);
    cl.try_parse(argc, argv);

    // check options
    if (!cl.exist("w-detector") && !cl.exist("d-detector") &&
        !cl.exist("n-detectors") && !cl.exist("small") && !cl.exist("big")) {
      throw(
          "need to specify either --w-detector, --d-detector or --n-detectors "
          "or --small or --big");
    }

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    if (cl.exist("small"))
      PET2D::Barrel::set_small_barrel_options(cl);
    else if (cl.exist("big"))
      PET2D::Barrel::set_big_barrel_options(cl);
    else
      PET2D::Barrel::calculate_scanner_options(cl);

    const auto& shape = cl.get<std::string>("shape");
    const auto& model_name = cl.get<std::string>("model");
    const auto& length_scale = cl.get<double>("base-length");

    Phantom phantom;
    // Read phantom
    for (auto& fn : cl.rest()) {
      std::ifstream in_phantom(fn);
      phantom.read_from_stream(in_phantom);
    }
    phantom.calculate_cdf();

    // run simmulation on given detector model & shape
    if (model_name == "always") {
      Common::AlwaysAccept<F> model;
      if (shape == "square") {
        run<SquareScanner>(cl, phantom, model);
      } else if (shape == "circle") {
        run<CircleScanner>(cl, phantom, model);
      } else if (shape == "triangle") {
        run<TriangleScanner>(cl, phantom, model);
      } else if (shape == "hexagon") {
        run<HexagonalScanner>(cl, phantom, model);
      }
    } else if (model_name == "scintillator") {
      Common::ScintillatorAccept<F> model(length_scale);
      if (shape == "square") {
        run<SquareScanner>(cl, phantom, model);
      } else if (shape == "circle") {
        run<CircleScanner>(cl, phantom, model);
      } else if (shape == "triangle") {
        run<TriangleScanner>(cl, phantom, model);
      } else if (shape == "hexagon") {
        run<HexagonalScanner>(cl, phantom, model);
      }
    }

    return 0;
  } catch (cmdline::exception& ex) {
    if (ex.help()) {
      std::cerr << ex.usage();
    }
    for (auto& msg : ex.errors()) {
      auto name = ex.name();
      if (name) {
        std::cerr << "error at " << name << ": " << msg << std::endl;
      } else {
        std::cerr << "error: " << msg << std::endl;
      }
    }
  } catch (std::string& ex) {
    std::cerr << "error: " << ex << std::endl;
  } catch (const char* ex) {
    std::cerr << "error: " << ex << std::endl;
  }
  return 1;
}

template <class DetectorClass, class PhantomClass, class ModelClass>
void run(cmdline::parser& cl, PhantomClass& phantom, ModelClass& model) {
  using MonteCarlo = Common::PhantomMonteCarlo<PhantomClass, DetectorClass>;
  using RNG = typename PhantomClass::RNG;

  auto& n_emissions = cl.get<int>("n-emissions");

  auto verbose = cl.exist("verbose");

  auto scanner = ScannerBuilder<DetectorClass>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, typename DetectorClass::F));
  scanner.set_sigma_dl(cl.get<float>("sigma"));
  if (cl.exist("tof-step"))
    scanner.set_tof_step(cl.get<double>("tof-step"));

  MonteCarlo monte_carlo(phantom, scanner);

  std::random_device rd;
  RNG rng(rd());
  if (cl.exist("seed")) {
    rng.seed(cl.get<std::mt19937::result_type>("seed"));
  }

  if (cl.exist("output")) {
    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    auto n_pixels = cl.get<int>("n-pixels");
    auto s_pixel = cl.get<double>("s-pixel");
    float ll = -s_pixel * n_pixels / 2;
    PET2D::PixelGrid<F, S> pixel_grid(
        n_pixels, n_pixels, s_pixel, PET2D::Point<F>(ll, ll));

    std::vector<int> image_emitted(n_pixels * n_pixels, 0);

    if (cl.exist("bin")) {
      int n_tof_positions = scanner.n_tof_positions(scanner.tof_step_size(),
                                                    scanner.max_dl_error());
      if (n_tof_positions == 0)
        n_tof_positions = 1;
      int n_detectors = scanner.size();
      std::vector<int> hits(n_detectors * n_detectors * n_tof_positions, 0);
      monte_carlo(
          rng,
          model,
          n_emissions,
          [&](const typename MonteCarlo::Event& event) {
            auto pixel = pixel_grid.pixel_at(event.center);
            auto pixel_index = pixel.y * n_pixels + pixel.x;
            image_emitted[pixel_index]++;
          },
          [&](const typename MonteCarlo::Event&,
              const typename MonteCarlo::FullResponse& full_response) {
            auto response = scanner.response_w_error(rng, full_response);
            if (response.tof_position < 0)
              response.tof_position = 0;
            if (response.tof_position >= n_tof_positions)
              response.tof_position = n_tof_positions - 1;
            int index = response.lor.first * n_detectors * n_tof_positions +
                        response.lor.second * n_tof_positions +
                        response.tof_position;
            hits[index]++;
          });

      std::ofstream out_hits(output);
      for (int d1 = 0; d1 < n_detectors; d1++)
        for (int d2 = 0; d2 < n_detectors; d2++)
          for (int tof = 0; tof < n_tof_positions; tof++) {
            if (hits[d1 * n_detectors * n_tof_positions + d2 * n_tof_positions +
                     tof] > 0)
              out_hits << d1 << " " << d2 << " " << tof << " "
                       << hits[d1 * n_detectors * n_tof_positions +
                               d2 * n_tof_positions + tof]
                       << "\n";
          }
    } else {
      std::ofstream out_wo_error(output_base_name + "_geom_only" + ext);
      std::ofstream out_w_error(output);
      std::ofstream out_exact_events(output_base_name + "_exact_events" + ext);
      std::ofstream out_full_response(output_base_name + "_full_response" +
                                      ext);

      monte_carlo(
          rng,
          model,
          n_emissions,
          [&](const typename MonteCarlo::Event& event) {
            auto pixel = pixel_grid.pixel_at(event.center);
            auto pixel_index = pixel.y * n_pixels + pixel.x;
            image_emitted[pixel_index]++;
          },
          [&](const typename MonteCarlo::Event& event,
              const typename MonteCarlo::FullResponse& full_response) {
            out_exact_events << event << "\n";
            out_full_response << full_response << "\n";
            out_wo_error << scanner.response_wo_error(full_response) << "\n";
            out_w_error << scanner.response_w_error(rng, full_response) << "\n";
          });
    }
    if (verbose) {
      std::cerr << "detected: " << monte_carlo.n_events_detected() << " events"
                << std::endl;
    }

    util::png_writer png(output_base_name + ".png");
    png.write(n_pixels, n_pixels, image_emitted);
  }
}

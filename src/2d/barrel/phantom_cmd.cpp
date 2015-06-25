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
#include "util/json_ostream.h"
#include "options.h"

#include "2d/strip/phantom.h"
#include "common/model.h"
#include "common/phantom_monte_carlo.h"

#if _OPENMP
#include <omp.h>
#endif

const double RADIAN = M_PI / 180;

using F = float;
using S = short;
using Hit = int;

using Pixel = PET2D::Pixel<S>;
using Point = PET2D::Point<F>;
using Event = PET2D::Event<F>;

template <typename DetectorType>
using Scanner = PET2D::Barrel::GenericScanner<DetectorType, MAX_DETECTORS, S>;
template <typename DetectorType>
using ScannerBuilder = PET2D::Barrel::ScannerBuilder<DetectorType>;

// all available detector shapes
using SquareScanner = Scanner<PET2D::Barrel::SquareDetector<F>>;
using CircleScanner = Scanner<PET2D::Barrel::CircleDetector<F>>;
using TriangleScanner = Scanner<PET2D::Barrel::TriangleDetector<F>>;
using HexagonalScanner = Scanner<PET2D::Barrel::PolygonalDetector<6, F>>;

using Ellipse = PET2D::Ellipse<F>;
using PhantomRegion = PET2D::Strip::PhantomRegion<F>;

template <typename DetectorType, typename Phantom, typename ModelType>
void run(cmdline::parser& cl, Phantom& phantom, ModelType& model);

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    PET2D::Barrel::add_phantom_options(cl);
    cl.add("small", 0, "small barrel", false);
    cl.add("big", 0, "big barrel", false);
    cl.add<float>("sigma", 0, "tof sigma", false, 0.06);
    cl.add("bin", 0, "ouput number of hits in each lor position");
    cl.try_parse(argc, argv);

#if _OPENMP
    if (cl.exist("n-threads")) {
      omp_set_num_threads(cl.get<int>("n-threads"));
    }
#endif

    if (cl.exist("big"))
      PET2D::Barrel::set_big_barrel_options(cl);
    PET2D::Barrel::calculate_scanner_options(cl);

    const auto& shape = cl.get<std::string>("shape");
    const auto& model_name = cl.get<std::string>("model");
    const auto& length_scale = cl.get<double>("base-length");

    auto verbose = cl.exist("verbose");

    std::vector<PhantomRegion> ellipse_list;

    // Read phantom
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

    PET2D::Strip::Phantom<F, S> phantom(ellipse_list);

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

template <typename DetectorType, typename Phantom, typename ModelType>
void run(cmdline::parser& cl, Phantom& phantom, ModelType& model) {

  auto& n_emissions = cl.get<int>("n-emissions");

  auto verbose = cl.exist("verbose");

  // NOTE: detector height will be determined per shape

  std::random_device rd;
  std::mt19937 gen(rd());
  if (cl.exist("seed")) {
    gen.seed(cl.get<std::mt19937::result_type>("seed"));
  }

  auto dr = ScannerBuilder<DetectorType>::build_multiple_rings(
      PET2D_BARREL_SCANNER_CL(cl, typename DetectorType::F));
  dr.set_sigma_dl(cl.get<float>("sigma"));
  if (cl.exist("tof-step"))
    dr.set_tof_step(cl.get<double>("tof-step"));

  Common::PhantomMonteCarlo<Phantom, DetectorType> monte_carlo(phantom, dr);

  typename Phantom::RNG rng;

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

  std::vector<int> hits;
  monte_carlo.generate(rng, model, n_emissions);
  if (cl.exist("bin")) {
    std::cerr << "bin\n";
    auto n_tof_positions =
        dr.n_tof_positions(dr.tof_step_size(), dr.max_dl_error());
    if (n_tof_positions == 0)
      n_tof_positions = 1;
    auto n_detectors = dr.size();
    hits.assign(n_detectors * n_detectors * n_tof_positions, 0);
    for (auto& full_response : monte_carlo) {
      auto response = dr.response_w_error(rng, full_response);
      if (response.tof_position < 0)
        response.tof_position = 0;
      if (response.tof_position >= n_tof_positions)
        response.tof_position = n_tof_positions - 1;
      int index = response.lor.first * n_detectors * n_tof_positions +
                  response.lor.second * n_tof_positions + response.tof_position;
      hits[index]++;
    }

    for (int d1 = 0; d1 < n_detectors; d1++)
      for (int d2 = 0; d2 < n_detectors; d2++)
        for (int tof = 0; tof < n_tof_positions; tof++) {
          if (hits[d1 * n_detectors * n_tof_positions + d2 * n_tof_positions +
                   tof] > 0)
            out_w_error << d1 << " " << d2 << " " << tof << " "
                        << hits[d1 * n_detectors * n_tof_positions +
                                d2 * n_tof_positions + tof] << "\n";
        }

  } else
    monte_carlo.write_out(rng);

  if (verbose) {
    std::cerr << "detected: " << monte_carlo.n_events_detected() << " events"
              << std::endl;
  }
}

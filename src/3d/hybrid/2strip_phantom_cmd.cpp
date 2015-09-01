#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/json.h"
#include "util/backtrace.h"
#include "util/progress.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "3d/geometry/phantom.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"
#include "3d/geometry/phantom_builder.h"

#include "scanner.h"

#include "common/model.h"
#include "common/phantom_monte_carlo.h"
#include "common/types.h"

using RNG = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, short, 8>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<RNG, F>;
using Allways = Common::AlwaysAccept<F>;
using Scintillator = Common::ScintillatorAccept<F>;
using Point = PET3D::Point<F>;
using Vector = PET3D::Vector<F>;
using Voxel = PET3D::Voxel<S>;
using Image = PET3D::VoxelMap<Voxel, F>;
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner, Image>;
using Event = MonteCarlo::Event;
using FullResponse = MonteCarlo::FullResponse;

namespace {
F strip_width = 0.005;
F strip_height = 0.019;
F strip_distance = 0.410;
F inner_radius = (strip_distance - strip_height) / 2;
F strip_length = 0.300;
}

int main(int argc, char* argv[]) {
  CMDLINE_TRY

  cmdline::parser cl;
  cl.footer("phantom_description.json ...");
  cl.add<cmdline::path>("output",
                        'o',
                        "output events file",
                        false,
                        "phantom.txt",
                        cmdline::not_from_file);
  cl.add<int>("n-emissions", 'e', "number of emmisions", false, 0);
  cl.add<double>("s-z", 0, "TOF sigma along z axis", false, 0.015);
  cl.add<double>("s-dl", 0, "TOF sigma delta-l", false, 0.06);
  cl.add<double>("radius", 'r', "cylinder radius", false, 0.0015);
  cl.add<double>("height", 'h', "cylinder height", false, 0.0020);
  cl.add<double>("x", 'x', "cylinder center x", false, 0);
  cl.add<double>("y", 'y', "cylinder center y", false, 0);
  cl.add<double>("z", 'z', "cylinder center z", false, 0);
  cl.add("verbose", 'v', "prints the iterations information on std::out");

  cl.parse_check(argc, argv);

  if (!cl.rest().size()) {
    if (argc == 1) {
      std::cerr << cl.usage();
      exit(0);
    } else {
      throw("at least one input phantom description expected, consult --help");
    }
  }

  auto verbose = cl.count("verbose");

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);
  scanner.set_sigmas(cl.get<double>("s-z"), cl.get<double>("s-dl"));

  using RNG = std::mt19937;
  RNG rng;
  Phantom::RegionPtrList regions;

  for (const auto& fn : cl.rest()) {
    std::ifstream in(fn);
    if (!in.is_open()) {
      throw("could not open file: " + cl.get<std::string>("phantoms"));
    }
    json j;
    j << in;

    if (!j.is_object()) {
      throw("no JSON object in file:" + cl.get<std::string>("phantoms"));
    }

    const json& j_phantoms = j["phantoms"];
    if (!j_phantoms.is_array()) {
      throw("phantoms array missing in JSON file: " +
            cl.get<std::string>("phantoms"));
    }

    for (const auto& j_phantom : j_phantoms) {
      auto region = PET3D::create_phantom_region_from_json<RNG, F>(j_phantom);
#if DEBUG
      std::cerr << "Adding region\n";
#endif
      regions.push_back(region);
    }
  }

  auto n_emissions = cl.get<int>("n-emissions");

  Phantom phantom(regions);

  Scintillator scintillator(0.100);
  MonteCarlo monte_carlo(phantom, scanner);

  if (cl.exist("output")) {
    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    std::ofstream out_json(output_base_name + ".json");
    out_json << json(scanner.barrel);

    std::ofstream out_wo_error(output_base_name + "_geom_only" + ext);
    std::ofstream out_w_error(output);
    std::ofstream out_exact_events(output_base_name + "_exact_events" + ext);
    std::ofstream out_full_response(output_base_name + "_full_response" + ext);

    util::progress progress(verbose, n_emissions, 10000);
    monte_carlo(
        rng,
        scintillator,
        n_emissions,
        [](const Event&) {},
        [&](const Event& event, const FullResponse& full_response) {
          out_exact_events << event << "\n";
          out_full_response << full_response << "\n";
          out_wo_error << scanner.response_wo_error(full_response) << "\n";
          out_w_error << scanner.response_w_error(rng, full_response) << "\n";
        },
        progress);
  }

  CMDLINE_CATCH
}

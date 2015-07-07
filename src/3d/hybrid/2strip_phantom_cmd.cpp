

#include "common/phantom_monte_carlo.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "common/model.h"

#include "scanner.h"

#include "3d/geometry/phantom.h"
#include "common/phantom_monte_carlo.h"

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/mathematica_ostream.h"
#include "util/json.h"
#include "util/json_ostream.h"

#include "3d/geometry/phantom_builder.h"

using F = float;
using S = short;
using RNG = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<F>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, short, 8>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<F, S, RNG>;
using Allways = Common::AlwaysAccept<F>;
using Scintillator = Common::ScintillatorAccept<F>;
using Point = PET3D::Point<F>;
using Vector = PET3D::Vector<F>;

using json = nlohmann::json;

namespace {
F strip_width = 0.005;
F strip_height = 0.019;
F strip_distance = 0.410;
F inner_radius = (strip_distance - strip_height) / 2;
F strip_length = 0.300;
}

int main(int argc, char* argv[]) {

  cmdline::parser cl;
  cl.add<cmdline::path>("output",
                        'o',
                        "output events file",
                        false,
                        "phantom.txt",
                        cmdline::not_from_file);
  cl.add<int>("n-emissions", 'e', "number of emmisions", false, 0);
  cl.add<float>("sigma-z", 0, "sigma-z", false, 0.015);
  cl.add<float>("sigma-dl", 0, "sigma-dl", false, 0.060);
  cl.add<float>("radius", 'r', "cylinder radius", false, 0.0015);
  cl.add<float>("height", 'h', "cylinder height", false, 0.0020);
  cl.add<float>("x", 'x', "cylinder center x", false, 0);
  cl.add<float>("y", 'y', "cylinder center y", false, 0);
  cl.add<float>("z", 'z', "cylinder center z", false, 0);
  cl.add<std::string>(
      "phantoms", '\0', "phantom description in JSON fromat", false);

  cl.parse_check(argc, argv);

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);
  scanner.set_sigmas(cl.get<float>("sigma-z"), cl.get<float>("sigma-dl"));

  util::mathematica_ostream mathematica(output_base_name + ".m");
  mathematica << scanner.barrel;

  util::json_ostream out_json(output_base_name + ".json");
  out_json << scanner.barrel;

  using RNG = std::mt19937;
  RNG rng;
  std::vector<PET3D::PhantomRegion<float, RNG>*> regions;

  if (cl.exist("phantoms")) {
    std::ifstream in(cl.get<std::string>("phantoms"));
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
      auto region =
          PET3D::create_phantom_region_from_json<float, RNG>(j_phantom);
#if DEBUG
      std::cerr << "Adding region\n";
#endif
      regions.push_back(region);
    }

  } else {
    float angle = std::atan2(0.0025f, 0.190);
    auto cylinder = new PET3D::CylinderRegion<float, RNG>(
        cl.get<float>("radius"),
        cl.get<float>("height"),
        1,
        PET3D::SphericalDistribution<float>(-angle, angle));
    PET3D::Matrix<float> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

    auto rotated_cylinder =
        new PET3D::RotatedPhantomRegion<float, RNG>(cylinder, R);
    Vector translation(
        cl.get<float>("x"), cl.get<float>("y"), cl.get<float>("z"));

    auto translated_cylinder = new PET3D::TranslatedPhantomRegion<float, RNG>(
        rotated_cylinder, translation);
    regions.push_back(translated_cylinder);
  }

  PET3D::Phantom<float, short, RNG> phantom(regions);

  Scintillator scintillator(0.100);
  Common::PhantomMonteCarlo<Phantom, Scanner> monte_carlo(phantom, scanner);

  std::ofstream out_wo_error(output_base_name + "_geom_only" + ext);
  monte_carlo.out_wo_error = out_wo_error;

  std::ofstream out_w_error(output);
  monte_carlo.out_w_error = out_w_error;

  std::ofstream out_exact_events(output_base_name + "_exact_events" + ext);
  monte_carlo.out_exact_events = out_exact_events;

  std::ofstream out_full_response(output_base_name + "_full_response" + ext);
  monte_carlo.out_full_response = out_full_response;
  monte_carlo.generate(rng, scintillator, cl.get<int>("n-emissions"));

  return 0;
}

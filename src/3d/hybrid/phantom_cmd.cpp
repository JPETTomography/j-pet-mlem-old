

#include "common/phantom_monte_carlo.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "common/model.h"

#include "scanner.h"

#include "3d/geometry/phantom.h"

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
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, S>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<RNG, F>;
using Allways = Common::AlwaysAccept<F>;
using Scintillator = Common::ScintillatorAccept<F>;
using Point = PET3D::Point<F>;
using Vector = PET3D::Vector<F>;

using json = nlohmann::json;

// FIXME: I don't know what is the purpose of this, but these are unused, so
// either should be removed or applied to the code.
#if HARDCODED_VALUES
namespace {
FType strip_width = 0.005;
FType strip_height = 0.019;
FType strip_distance = 0.410;
FType inner_radius = (strip_distance - strip_height) / 2;
FType strip_length = 0.300;
}
#endif

int main(int argc, char* argv[]) {

  try {
    cmdline::parser cl;
    // cl.add<int>("n-emissions", 'e', "number of emmisions", false, 0);
    cl.add<float>("sigma-z", 0, "sigma-z", false, 0.015);
    cl.add<float>("sigma-dl", 0, "sigma-dl", false, 0.060);
    cl.add("small", 0, "small barrel");
    cl.add("big", 0, "big barrel");
    cl.add<std::string>(
        "phantoms", '\0', "phantom description in JSON format", true);
    cl.add<double>("z-position", 'z', "position of the z plane", false, 0);
    cl.add<double>("length", 0, "length of the detector", false, 0.3);

    PET3D::Hybrid::add_phantom_options(cl);

    cl.parse_check(argc, argv);

    // check options
    if (!cl.exist("w-detector") && !cl.exist("d-detector") &&
        !cl.exist("n-detectors") && !cl.exist("small") && !cl.exist("big")) {
      throw(
          "need to specify either --w-detector, --d-detector or --n-detectors "
          "or --small or --big");
    }

    if (cl.exist("small"))
      PET3D::Hybrid::set_small_barrel_options(cl);
    else if (cl.exist("big"))
      PET3D::Hybrid::set_big_barrel_options(cl);
    else
      PET3D::Hybrid::calculate_scanner_options(cl);

    auto output = cl.get<cmdline::path>("output");
    auto output_base_name = output.wo_ext();
    auto ext = output.ext();

    Scanner scanner = Scanner::build_scanner_from_cl(cl);
    scanner.set_sigmas(cl.get<float>("sigma-z"), cl.get<float>("sigma-dl"));

    util::mathematica_ostream mathematica(output_base_name + ".m");
    mathematica << scanner.barrel;

    util::json_ostream out_json(output_base_name + ".json");
    out_json << scanner.barrel;

    using RNG = std::mt19937;
    RNG rng;
    Phantom::RegionPtrList regions;

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
        auto region = PET3D::create_phantom_region_from_json<RNG, F>(j_phantom);
#if DEBUG
        std::cerr << "Adding region\n";
#endif
        regions.push_back(region);
      }

    } else {
      float angle = std::atan2(0.0025f, 0.190);
      auto cylinder = new Phantom::CylinderRegion<>(
          cl.get<float>("radius"),
          cl.get<float>("height"),
          1,
          PET3D::Distribution::SphericalDistribution<float>(-angle, angle));
      PET3D::Matrix<float> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

      auto rotated_cylinder = new Phantom::RotatedRegion(cylinder, R);
      Vector translation(
          cl.get<float>("x"), cl.get<float>("y"), cl.get<float>("z"));

      auto translated_cylinder =
          new Phantom::TranslatedRegion(rotated_cylinder, translation);
      regions.push_back(translated_cylinder);
    }

    Phantom phantom(regions);

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
    monte_carlo.write_out(rng);

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

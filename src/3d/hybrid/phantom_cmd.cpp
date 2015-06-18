

#include "phantom_monte_carlo.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/model.h"

#include "scanner.h"

#include "3d/geometry/phantom.h"
#include "phantom_monte_carlo.h"

#include "cmdline.h"
#include "util/cmdline_types.h"
#include "util/cmdline_hooks.h"
#include "util/mathematica_ostream.h"
#include "util/json_ostream.h"

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include "3d/geometry/phantom_builder.h"

using FType = float;
using SType = short;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 192, short>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<FType, SType, RNGType>;
using Allways = PET2D::Barrel::AlwaysAccept<FType>;
using Scintillator = PET2D::Barrel::ScintillatorAccept<FType>;
using Point = PET3D::Point<FType>;
using Vector = PET3D::Vector<FType>;

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

  if (!cl.exist("small") && !cl.exist("big")) {
    throw("need to specify either ---small or --big");
  }

  if (cl.exist("small"))
    PET3D::Hybrid::set_small_barrel_options(cl);
  if (cl.exist("big"))
    PET3D::Hybrid::set_big_barrel_options(cl);

  PET3D::Hybrid::calculate_scanner_options(cl);
  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();

  Scanner scanner = Scanner::build_scanner_from_cl(cl);
  scanner.set_sigmas(cl.get<float>("sigma-z"), cl.get<float>("sigma-dl"));

  util::mathematica_ostream mathematica(output_base_name + ".m");
  mathematica << scanner.barrel;

  util::json_ostream json(output_base_name + ".json");
  json << scanner.barrel;

  using RNG = std::mt19937;
  RNG rng;
  std::vector<PET3D::PhantomRegion<float, RNG>*> regions;

  if (cl.exist("phantoms")) {
    FILE* in = fopen(cl.get<std::string>("phantoms").c_str(), "r");
    if (!in) {
      std::cerr << "could not open file: `" << cl.get<std::string>("phantoms")
                << "'\n";
      exit(0);
    }
    char readBuffer[256];
    rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
    rapidjson::Document doc;
    doc.ParseStream(input_stream);

    if (!doc.IsObject()) {
      std::cerr << "file `" << cl.get<std::string>("phantoms")
                << "' does not contain a JSON object " << doc.GetType() << "\n";
      exit(0);
    }
    const rapidjson::Value& phantoms_val = doc["phantoms"];
    if (!phantoms_val.IsArray()) {
      std::cerr << "file `" << cl.get<std::string>("phantoms")
                << "' does not contain Phantoms\n";
      exit(0);
    }

    for (auto it = phantoms_val.Begin(); it != phantoms_val.End(); it++) {

      auto region = PET3D::create_phantom_region_from_json<float, RNG>(*it);
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

  Allways allways;
  Scintillator scintillator(0.100);
  PET3D::PhantomMonteCarlo<Phantom, Scanner> monte_carlo(phantom, scanner);

  std::ofstream no_error_stream(output_base_name + "_geom_only" + ext);
  monte_carlo.set_no_error_stream(no_error_stream);

  std::ofstream error_stream(output);
  monte_carlo.set_error_stream(error_stream);

  std::ofstream exact_event_stream(output_base_name + "_exact_events" + ext);
  monte_carlo.set_exact_event_stream(exact_event_stream);

  std::ofstream full_response_stream(output_base_name + "_full_response" + ext);
  monte_carlo.set_full_response_stream(full_response_stream);
  monte_carlo.generate(rng, scintillator, cl.get<int>("n-emissions"));

  return 0;
}



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

using FType = float;
using SType = short;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 8, short>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<FType, SType, RNGType>;
using Allways = PET2D::Barrel::AlwaysAccept<FType>;
using Scintillator = PET2D::Barrel::ScintillatorAccept<FType>;
using Point = PET3D::Point<FType>;
using Vector = PET3D::Vector<FType>;

namespace {
FType strip_width = 0.005;
FType strip_height = 0.019;
FType strip_distance = 0.410;
FType inner_radius = (strip_distance - strip_height) / 2;
FType strip_length = 0.300;
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
  cl.parse_check(argc, argv);

  auto output = cl.get<cmdline::path>("output");
  auto output_base_name = output.wo_ext();
  auto ext = output.ext();

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);
  scanner.set_sigmas(cl.get<float>("sigma-z"), cl.get<float>("sigma-dl"));

  using RNG = std::mt19937;
  RNG rng;
  std::vector<PET3D::PhantomRegion<float, RNG>*> regions;
  float angle = std::atan2(0.0025f, 0.190);
  auto cylinder = new PET3D::CylinderRegion<float, RNG>(
      0.0025, 0.001, 1, PET3D::SphericalDistribution<float>(-angle, angle));
  PET3D::Matrix<float> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

  auto rotated_cylinder =
      new PET3D::RotatedPhantomRegion<float, RNG>(cylinder, R);
  regions.push_back(rotated_cylinder);
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

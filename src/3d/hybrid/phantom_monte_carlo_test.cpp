#include "util/test.h"

#include "phantom_monte_carlo.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "2d/barrel/model.h"

#include "scanner.h"

#include "3d/geometry/phantom.h"
#include "phantom_monte_carlo.h"

using FType = float;
using SType = short;
using RNGType = std::mt19937;
using Detector = PET2D::Barrel::SquareDetector<FType>;
using Scanner2D = PET2D::Barrel::GenericScanner<Detector, 8, short>;
using Scanner = PET3D::Hybrid::Scanner<Scanner2D>;
using Phantom = PET3D::Phantom<FType, SType, RNGType>;
using Allways = PET2D::Barrel::AlwaysAccept<FType>;
using Scintillator = PET2D::Barrel::ScintillatorAccept<FType>;
namespace {
FType strip_width = 0.005;
FType strip_height = 0.019;
FType strip_distance = 0.410;
FType inner_radius = (strip_distance - strip_height) / 2;
FType strip_length = 0.300;
}

TEST("PET3D/hubrid/phantom_monte_carlo") {

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);

  using RNG = std::mt19937;
  RNG rng;
  std::vector<PET3D::PhantomRegion<float, RNG>*> regions;
  float angle = std::atan2(0.0025f, 0.400f);
  auto cylinder = new PET3D::CylinderRegion<float, RNG>(
      0.0015, 0.001, 1, PET3D::SphericalDistribution<float>(-angle, angle));
  PET3D::Matrix<float> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

  auto rotated_cylinder =
      new PET3D::RotatedPhantomRegion<float, RNG>(cylinder, R);
  regions.push_back(rotated_cylinder);
  PET3D::Phantom<float, short, RNG> phantom(regions);

  Allways allways;

  PhantomMonteCarlo<Phantom, Scanner> monte_carlo(phantom, scanner);
  monte_carlo.generate(rng, allways, 1000);
}
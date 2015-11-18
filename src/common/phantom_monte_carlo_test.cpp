#include "util/test.h"

#include "2d/barrel/square_detector.h"
#include "2d/barrel/generic_scanner.h"
#include "2d/barrel/scanner_builder.h"
#include "3d/hybrid/scanner.h"
#include "3d/geometry/phantom.h"
#include "3d/geometry/voxel.h"
#include "3d/geometry/voxel_map.h"

#include "phantom_monte_carlo.h"

#include "common/model.h"
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
using MonteCarlo = Common::PhantomMonteCarlo<Phantom, Scanner>;

namespace {
F strip_width = 0.005;
F strip_height = 0.019;
F strip_distance = 0.410;
F inner_radius = (strip_distance - strip_height) / 2;
F strip_length = 0.300;
}

TEST("common/phantom_monte_carlo/point_source") {

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);
  scanner.set_sigmas(0.010, 0.024);

  using RNG = std::mt19937;
  RNG rng;
  Phantom::RegionPtrList regions;

  auto emitter = new Phantom::PointRegion<
      PET3D::Distribution::SingleDirectionDistribution<F>>(
      1,
      PET3D::Distribution::SingleDirectionDistribution<F>(
          Vector::from_euler_angles(0, 2 * M_PI / 6)),
      Point(0, 0, 0));

  regions.push_back(emitter);
  Phantom phantom(regions);

  Scintillator scintillator(0.100);
  MonteCarlo monte_carlo(phantom, scanner);

  monte_carlo(rng,
              scintillator,
              10000,
              [](Phantom::Event&) {},
              [](Phantom::Event&, Scanner::FullResponse&) {},
              [](int, bool) {});

  // FIXME: put there some better estimation how many events we shall catch
  REQUIRE(monte_carlo.n_events_detected() > 100);
}

TEST("common/phantom_monte_carlo/phantom_region") {

  auto scanner2d = PET2D::Barrel::ScannerBuilder<Scanner2D>::build_single_ring(
      inner_radius, 2, strip_width, strip_height);

  Scanner scanner(scanner2d, strip_length);
  scanner.set_sigmas(0.010, 0.024);

  using RNG = std::mt19937;
  RNG rng;
  Phantom::RegionPtrList regions;
  F angle = std::atan2(0.0025, 0.400);
  auto cylinder = new Phantom::CylinderRegion<>(
      0.0015,
      0.001,
      1,
      PET3D::Distribution::SphericalDistribution<F>(-angle, angle));
  PET3D::Matrix<F> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

  auto rotated_cylinder = new Phantom::RotatedRegion(cylinder, R);
  regions.push_back(rotated_cylinder);
  Phantom phantom(regions);

  Scintillator scintillator(0.100);
  MonteCarlo monte_carlo(phantom, scanner);

  monte_carlo(rng,
              scintillator,
              10000,
              [](Phantom::Event&) {},
              [](Phantom::Event&, Scanner::FullResponse&) {},
              [](int, bool) {});

  // FIXME: put there some better estimation how many events we shall catch
  REQUIRE(monte_carlo.n_events_detected() > 100);
}

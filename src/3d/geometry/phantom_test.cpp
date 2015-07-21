#include <random>
#include <iostream>
#include <fstream>

#include "util/test.h"

#include "matrix.h"
#include "phantom.h"

#include "common/types.h"

using RNG = std::mt19937;
using Phantom = PET3D::Phantom<RNG, F>;
using Point = Phantom::Point;
using AngularDistribution = PET3D::Distribution::SphericalDistribution<F>;

TEST("3d/geometry/phantom/cylinder_region") {
  using Region = Phantom::CylinderRegion<AngularDistribution>;

  Region region(2, 3, 1, AngularDistribution(-M_PI / 3, M_PI / 3));

  REQUIRE(region.volume() == Approx(4 * M_PI * 3).epsilon(1e-7));
  REQUIRE(region.intensity == 1.0_e7);

  Point p1(1.2, 0.1, 1.4);
  REQUIRE(region.in(p1));
  Point p2(1.2, 0.0, 1.7);
  REQUIRE(!region.in(p2));
  Point p3(-2.1, 0.05, -1.0);
  REQUIRE(!region.in(p3));

  std::mt19937 rng;

  for (int i = 0; i < 100; i++) {
    auto event = region.random_event(rng);
    auto direction = event.direction;

    REQUIRE(((direction.z <= std::sqrt(3) / 2) &&
             (direction.z >= -std::sqrt(3) / 2)));
  }
}

TEST("3d/geometry/phantom/cylinder") {
  RNG rng;
  Phantom::RegionPtrList regions;
  F angle = std::atan2(0.0025, 0.400);
  auto cylinder = new Phantom::CylinderRegion<AngularDistribution>(
      0.0015, 0.001, 1, AngularDistribution(-angle, angle));
  PET3D::Matrix<F> R{ 1, 0, 0, 0, 0, 1, 0, 1, 0 };

  auto rotated_cylinder = new Phantom::RotatedRegion(cylinder, R);
  regions.push_back(rotated_cylinder);
  Phantom phantom(regions);

  for (int i = 0; i < 10000; i++) {
    auto event = phantom.gen_event(rng);
    auto p = event.origin;
    (void)p;  // FIXME: test position here
  }
}

TEST("3d/geometry/phantom/ellipsoid") {
  using RNG = std::mt19937;
  RNG rng;
  Phantom::RegionPtrList regions;
  F angle = std::atan2(0.0025, 0.400);
  auto ellipsoid = new Phantom::EllipsoidRegion<AngularDistribution>(
      0.005, 0.01, 0.02, 1, AngularDistribution(-angle, angle));

  regions.push_back(ellipsoid);
  Phantom phantom(regions);

  for (int i = 0; i < 10000; i++) {
    auto event = phantom.gen_event(rng);
    auto p = event.origin;
    (void)p;  // FIXME: test position here
  }
}

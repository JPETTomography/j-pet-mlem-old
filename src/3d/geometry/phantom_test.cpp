#include <random>
#include <iostream>
#include <fstream>

#include "util/test.h"

#include "matrix.h"
#include "phantom.h"

TEST("3d/geometry/phantom/cylinder_region") {
  using F = float;
  using Region =
      PET3D::CylinderRegion<F, std::mt19937, PET3D::SphericalDistribution<F>>;
  using Point = Region::Point;

  Region region(
      2.0f, 3.0f, 1.0f, PET3D::SphericalDistribution<F>(-M_PI / 3, M_PI / 3));

  REQUIRE(region.volume() == Approx(4.0 * M_PI * 3.0).epsilon(1e-7));
  REQUIRE(region.intensity == 1.0_e7);

  Point p1(1.2f, 0.1f, 1.4f);
  REQUIRE(region.in(p1));
  Point p2(1.2f, 0.0f, 1.7f);
  REQUIRE(!region.in(p2));
  Point p3(-2.1f, 0.05f, -1.0f);
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

  std::ofstream out("test_output/cylinder.txt");

  for (int i = 0; i < 10000; i++) {
    auto event = phantom.gen_event(rng);
    auto p = event.origin;
    auto v = event.direction;
    out << p.x << " " << p.y << " " << p.z << " ";
    out << v.x << " " << v.y << " " << v.z << "\n";
  }
  out.close();
}

TEST("3d/geometry/phantom/ellipsoid") {
  using RNG = std::mt19937;
  RNG rng;
  std::vector<PET3D::PhantomRegion<float, RNG>*> regions;
  float angle = std::atan2(0.0025f, 0.400f);
  auto ellipsoid = new PET3D::EllipsoidRegion<float, RNG>(
      0.005, 0.01, 0.02, 1, PET3D::SphericalDistribution<float>(-angle, angle));

  regions.push_back(ellipsoid);
  PET3D::Phantom<float, short, RNG> phantom(regions);

  std::ofstream out("test_output/ellipsoid.txt");

  for (int i = 0; i < 10000; i++) {
    auto event = phantom.gen_event(rng);
    auto p = event.origin;
    auto v = event.direction;
    out << p.x << " " << p.y << " " << p.z << " ";
    out << v.x << " " << v.y << " " << v.z << "\n";
  }
  out.close();
}

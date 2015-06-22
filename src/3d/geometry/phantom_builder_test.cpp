#include <iostream>
#include <fstream>

#include "util/test.h"
#include "util/json.h"

#include "3d/geometry/phantom_builder.h"
#include "3d/geometry/event_generator.h"

TEST("3d/geometry/phantom_builder/rapid_json") {
  std::ifstream in("src/3d/hybrid/point_source.json");
  if (!in.is_open()) {
    FAIL("cannot open src/3d/hybrid/point_source.json");
  }
  json j;
  j << in;

  REQUIRE(j.is_array());
  const json& j_obj = j[0];

  REQUIRE(j_obj.count("type"));
  REQUIRE(j_obj.count("angular"));
}

TEST("3d/geometry/phantom_builder/angular_distribution") {
  std::ifstream in("src/3d/hybrid/point_source.json");
  if (!in.is_open()) {
    FAIL("cannot open src/3d/hybrid/point_source.json");
  }
  json j;
  j << in;

  const json& j_obj = j[0];
  const json& j_angular = j_obj["angular"];

  PET3D::SingleDirectionDistribution<float> distr =
      PET3D::create_angular_distribution_from_json<
          PET3D::SingleDirectionDistribution<float>>(j_angular);

  REQUIRE(distr.direction.x == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
  REQUIRE(distr.direction.y == Approx(0.0f).epsilon(1e-7));
  REQUIRE(distr.direction.z == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));

  int dummy;
  auto dir = distr(dummy);

  REQUIRE(dir.x == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
  REQUIRE(dir.y == Approx(0.0f).epsilon(1e-7));
  REQUIRE(dir.z == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
}

TEST("3d/geometry/phantom_builder/angular_distribution/spherical",
     "spherical") {
  std::ifstream in("src/3d/geometry/test_phantoms.json");
  if (!in.is_open()) {
    FAIL("cannot open src/3d/geometry/test_phantoms.json");
  }
  json j;
  j << in;

  const json& j_phantoms = j["phantoms"];
  const json& j_phantom = j_phantoms[0];
  const json& j_angular = j_phantom["angular"];

  PET3D::SphericalDistribution<float> distr =
      PET3D::create_angular_distribution_from_json<
          PET3D::SphericalDistribution<float>>(j_angular);

  REQUIRE(distr.theta_min == -0.01_e7);
  REQUIRE(distr.theta_max == 0.01_e7);
}

TEST("3d/geometry/phantom_builder/phantom") {
  using RNGType = std::mt19937;
  std::ifstream in("src/3d/geometry/test_phantoms.json");
  if (!in.is_open()) {
    FAIL("cannot open src/3d/geometry/test_phantoms.json");
  }
  json j;
  j << in;

  const json& j_phantoms = j["phantoms"];
  REQUIRE(j_phantoms.is_array());

  {
    const json& j_phantom = j_phantoms[0];
    auto phantom = static_cast<PET3D::CylinderRegion<float, RNGType>*>(
        PET3D::create_phantom_region_from_json<float, RNGType>(j_phantom));

    REQUIRE(phantom->intensity == 1.0_e7);
    REQUIRE(phantom->radius == 0.005_e7);
    REQUIRE(phantom->height == 0.002_e7);
  }

  {
    const json& j_phantom = j_phantoms[1];
    auto phantom =
        PET3D::create_phantom_region_from_json<float, RNGType>(j_phantom);

    REQUIRE(phantom->in(PET3D::Point<float>(-0.05, 0.007, 0.03)));
    REQUIRE(!phantom->in(PET3D::Point<float>(-0.05, 0.011, 0.03)));
  }

  {
    const json& j_phantom = j_phantoms[2];
    auto phantom =
        PET3D::create_phantom_region_from_json<float, RNGType>(j_phantom);

    REQUIRE(phantom->in(PET3D::Point<float>(0.05, 0.0, -0.10)));
    REQUIRE(!phantom->in(PET3D::Point<float>(0.0, .16, 0.0)));
  }
}

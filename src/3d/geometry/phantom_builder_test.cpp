#include "util/test.h"
#include "util/json.h"

#include "3d/geometry/phantom_builder.h"
#include "3d/geometry/event_generator.h"

static const char* point_source_json = R"JSON([
  {
    "angular": {
      "direction": [
        1,
        0,
        1
      ],
      "type": "single-direction"
    },
    "intensity": 1.0,
    "origin": [0, 0, 0],
    "type": "point-source"
  }
])JSON";

TEST("3d/geometry/phantom_builder/rapid_json") {
  json j = json::parse(point_source_json);

  REQUIRE(j.is_array());
  const json& j_obj = j[0];

  REQUIRE(j_obj.count("type"));
  REQUIRE(j_obj.count("angular"));
}

TEST("3d/geometry/phantom_builder/angular_distribution") {
  json j = json::parse(point_source_json);

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

static const char* test_phantoms_json = R"JSON({
  "phantoms": [
    {
      "angular": {
        "theta-max": 0.01,
        "theta-min": -0.01,
        "type": "spherical"
      },
      "height": 0.002,
      "id": "cylinder",
      "intensity": 1.0,
      "radius": 0.005,
      "type": "cylinder"
    },
    {
      "displacement": [
        -0.05,
        0.0,
        0.03
      ],
      "phantom": {
        "R": [1, 0, 0,
              0, 0, 1,
              0, 1, 0],
        "phantom": {
          "angular": {
            "type": "spherical"
          },
          "height": 0.02,
          "id": "cylinder",
          "intensity": 1.0,
          "radius": 0.005,
          "type": "cylinder"
        },
        "type": "rotated"
      },
      "type": "translated"
    },
    {
      "angular": {
        "type": "spherical"
      },
      "rx": 0.1,
      "ry": 0.15,
      "rz": 0.2,
      "type": "ellipsoid"
    }
  ]
})JSON";

TEST("3d/geometry/phantom_builder/angular_distribution/spherical") {
  json j = json::parse(test_phantoms_json);

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
  using RNG = std::mt19937;
  using Phantom = PET3D::Phantom<RNG, float>;
  using Point = Phantom::Point;

  json j = json::parse(test_phantoms_json);

  const json& j_phantoms = j["phantoms"];
  REQUIRE(j_phantoms.is_array());

  {
    const json& j_phantom = j_phantoms[0];
    auto phantom = static_cast<Phantom::CylinderRegion<>*>(
        PET3D::create_phantom_region_from_json<RNG, float>(j_phantom));

    REQUIRE(phantom->intensity == 1.0_e7);
    REQUIRE(phantom->radius == 0.005_e7);
    REQUIRE(phantom->height == 0.002_e7);
  }

  {
    const json& j_phantom = j_phantoms[1];
    auto phantom =
        PET3D::create_phantom_region_from_json<RNG, float>(j_phantom);

    REQUIRE(phantom->in(Point(-0.05, 0.007, 0.03)));
    REQUIRE(!phantom->in(Point(-0.05, 0.011, 0.03)));
  }

  {
    const json& j_phantom = j_phantoms[2];
    auto phantom =
        PET3D::create_phantom_region_from_json<RNG, float>(j_phantom);

    REQUIRE(phantom->in(Point(0.05, 0.0, -0.10)));
    REQUIRE(!phantom->in(Point(0.0, .16, 0.0)));
  }
}

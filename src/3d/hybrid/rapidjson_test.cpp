#include <iostream>
#include <cstdio>

#include "util/test.h"
#include "3d/geometry/phantom_builder.h"

#include "rapidjson/document.h"
#include "rapidjson/filereadstream.h"

#include "3d/geometry/event_generator.h"
#include "3d/geometry/phantom_builder.h"

TEST("rapid_json") {
  FILE* in = fopen("src/3d/hybrid/point_source.js", "r");
  if (!in) {
    std::cerr << "cannot open src/3d/hybrid/point_source.js\n";
  }

  rapidjson::Document doc;
  char readBuffer[256];
  rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
  doc.ParseStream(input_stream);

  REQUIRE(doc.IsArray());
  rapidjson::Value& obj = doc[0];

  REQUIRE(obj.HasMember("type"));
  REQUIRE(obj.HasMember("angular"));
  fclose(in);
}

TEST("PhantomBuilder/angular_distribution") {
  FILE* in = fopen("src/3d/hybrid/point_source.js", "r");
  if (!in) {
    std::cerr << "cannot open src/3d/hybrid/point_source.js\n";
  }

  rapidjson::Document doc;
  char readBuffer[256];
  rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
  doc.ParseStream(input_stream);
  fclose(in);

  rapidjson::Value& obj = doc[0];
  rapidjson::Value& angular_val = obj["angular"];

  PET3D::SingleDirectionDistribution<float> distr =
      PET3D::create_angular_distribution_from_json<
          PET3D::SingleDirectionDistribution<float>>(angular_val);

  REQUIRE(distr.direction.x == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
  REQUIRE(distr.direction.y == Approx(0.0f).epsilon(1e-7));
  REQUIRE(distr.direction.z == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));

  int dummy;
  auto dir = distr(dummy);

  REQUIRE(dir.x == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
  REQUIRE(dir.y == Approx(0.0f).epsilon(1e-7));
  REQUIRE(dir.z == Approx(1.0f / std::sqrt(2.0f)).epsilon(1e-7));
}

TEST("PhantomBuilder/angular_distribution/spherical", "spherical") {
  FILE* in = fopen("src/3d/hybrid/disk.json", "r");
  if (!in) {
    std::cerr << "cannot open src/3d/hybrid/disk.json\n";
  }

  rapidjson::Document doc;
  char readBuffer[256];
  rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
  doc.ParseStream(input_stream);
  fclose(in);

  rapidjson::Value& phantoms = doc["phantoms"];
  rapidjson::Value& obj = phantoms[0];
  rapidjson::Value& angular_val = obj["angular"];

  PET3D::SphericalDistribution<float> distr =
      PET3D::create_angular_distribution_from_json<
          PET3D::SphericalDistribution<float>>(angular_val);

  REQUIRE(distr.theta_min == -0.01_e7);
  REQUIRE(distr.theta_max == 0.01_e7);
}

TEST("PhantomBuilder/phantom") {
  using RNGType = std::mt19937;
  FILE* in = fopen("src/3d/hybrid/disk.json", "r");
  if (!in) {
    std::cerr << "cannot open src/3d/hybrid/disk.json\n";
  }

  rapidjson::Document doc;
  char readBuffer[256];
  rapidjson::FileReadStream input_stream(in, readBuffer, sizeof(readBuffer));
  doc.ParseStream(input_stream);
  fclose(in);

  rapidjson::Value& phantoms_val = doc["phantoms"];
  {
    rapidjson::Value& phantom_val = phantoms_val[0];
    auto phantom = static_cast<PET3D::CylinderRegion<float, RNGType>*>(
        PET3D::create_phantom_region_from_json<float, RNGType>(phantom_val));

    REQUIRE(phantom->intensity == 1.0_e7);
    REQUIRE(phantom->radius == 0.005_e7);
    REQUIRE(phantom->height == 0.002_e7);
  }

  {
    rapidjson::Value& phantom_val = phantoms_val[1];
    auto phantom =
        PET3D::create_phantom_region_from_json<float, RNGType>(phantom_val);

    REQUIRE(phantom->in(PET3D::Point<float>(-0.05, 0.007, 0.03)));
    REQUIRE(!phantom->in(PET3D::Point<float>(-0.05, 0.011, 0.03)));
  }
}

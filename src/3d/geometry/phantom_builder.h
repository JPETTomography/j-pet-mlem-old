#pragma once

#include "point.h"
#include "vector"
#include "event_generator.h"
#include "phantom.h"

#include "util/json.h"

using json = nlohmann::json;

namespace PET3D {

template <typename AngularDistribution>
AngularDistribution create_angular_distribution_from_json(const json& obj);

template <>
SingleDirectionDistribution<float> create_angular_distribution_from_json<
    SingleDirectionDistribution<float>>(const json& j) {
  using Vector = Vector<float>;

  std::string type = j["type"];

  if (type != "single-direction") {
    throw("type mismatch in SingleDirectionDistribution: " + type);
  }

  const json& direction_val = j["direction"];
  Vector direction(direction_val[0], direction_val[1], direction_val[2]);
  return SingleDirectionDistribution<float>(direction);
}

template <>
SphericalDistribution<float> create_angular_distribution_from_json<
    SphericalDistribution<float>>(const json& j) {

  // std::cerr<<"creating  angular distribution\n";
  std::string type = j["type"];

  if (type != "spherical") {
    throw("type mismatch in SphericalDistribution: " + type);
  }

  float theta_min = -M_PI / 2;
  float theta_max = M_PI / 2;

  if (j.count("theta-min")) {
    theta_min = j["theta-min"];
  }
  if (j.count("theta-max")) {
    theta_max = j["theta-max"];
  }

  return SphericalDistribution<float>(theta_min, theta_max);
}

template <typename FType, typename RNG>
PhantomRegion<FType, RNG>* create_phantom_region_from_json(const json& j) {

  if (!j.count("type")) {
    throw("phantom region does not have type member");
  }

  std::string type = j["type"];

  if (type == "cylinder") {
    FType radius = j["radius"];
    FType height = j["height"];
    FType intensity = j["intensity"];

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      SphericalDistribution<FType> angular =
          create_angular_distribution_from_json<SphericalDistribution<FType>>(
              j_angular);
      return new CylinderRegion<FType, RNG>(radius, height, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "ellipsoid") {
    FType rx = j["rx"];
    FType ry = j["ry"];
    FType rz = j["rz"];
    FType intensity =
        j.count("intensity") ? j["intensity"].get<FType>() : FType(1);

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      SphericalDistribution<FType> angular =
          create_angular_distribution_from_json<SphericalDistribution<FType>>(
              j_angular);
      return new EllipsoidRegion<FType, RNG>(rx, ry, rz, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "point") {
    throw("not implemented yet");

  } else if (type == "rotated") {
    auto phantom = create_phantom_region_from_json<FType, RNG>(j["phantom"]);
    const json& j_R = j["R"];
    PET3D::Matrix<FType> R;
    int i = 0;
    for (const auto& el : j_R) {
      R(i++) = el;
    }
    return new RotatedPhantomRegion<FType, RNG>(phantom, R);

  } else if (type == "translated") {
    auto phantom = create_phantom_region_from_json<FType, RNG>(j["phantom"]);
    const json& displacement_val = j["displacement"];
    PET3D::Vector<FType> displacement;
    displacement.x = displacement_val[0];
    displacement.y = displacement_val[1];
    displacement.z = displacement_val[2];
    return new TranslatedPhantomRegion<FType, RNG>(phantom, displacement);

  } else {
    throw("unknown region type: " + type);
    return nullptr;
  }
}
}

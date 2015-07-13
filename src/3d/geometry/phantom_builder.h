#pragma once

#include "point.h"
#include "vector"
#include "event_generator.h"
#include "phantom.h"

#include "util/json.h"

namespace PET3D {

template <typename AngularDistribution>
AngularDistribution create_angular_distribution_from_json(const json& obj);

template <>
Distribution::SingleDirectionDistribution<float>
create_angular_distribution_from_json<
    Distribution::SingleDirectionDistribution<float>>(const json& j) {
  using Vector = Vector<float>;

  std::string type = j["type"];

  if (type != "single-direction") {
    throw("type mismatch in SingleDirectionDistribution: " + type);
  }

  const json& direction_val = j["direction"];
  Vector direction(direction_val[0], direction_val[1], direction_val[2]);
  return Distribution::SingleDirectionDistribution<float>(direction);
}

template <>
Distribution::SphericalDistribution<float>
create_angular_distribution_from_json<
    Distribution::SphericalDistribution<float>>(const json& j) {

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

  return Distribution::SphericalDistribution<float>(theta_min, theta_max);
}

template <class RNG, typename FType>
typename Phantom<RNG, FType>::Region* create_phantom_region_from_json(
    const json& j) {
  using F = FType;
  using AngularDistribution = Distribution::SphericalDistribution<F>;
  using Phantom = PET3D::Phantom<RNG, F>;
  using CylinderRegion =
      typename Phantom::template CylinderRegion<AngularDistribution>;
  using EllipsoidRegion =
      typename Phantom::template EllipsoidRegion<AngularDistribution>;
  using RotatedRegion = typename Phantom::RotatedRegion;
  using TranslatedRegion = typename Phantom::TranslatedRegion;

  if (!j.count("type")) {
    throw("phantom region does not have type member");
  }

  std::string type = j["type"];

  if (type == "cylinder") {
    F radius = j["radius"];
    F height = j["height"];
    F intensity = j["intensity"];

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      AngularDistribution angular =
          create_angular_distribution_from_json<AngularDistribution>(j_angular);
      return new CylinderRegion(radius, height, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "ellipsoid") {
    F rx = j["rx"];
    F ry = j["ry"];
    F rz = j["rz"];
    F intensity = j.count("intensity") ? j["intensity"].get<F>() : F(1);

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      AngularDistribution angular =
          create_angular_distribution_from_json<AngularDistribution>(j_angular);
      return new EllipsoidRegion(rx, ry, rz, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "point") {
    throw("not implemented yet");

  } else if (type == "rotated") {
    auto phantom = create_phantom_region_from_json<RNG, F>(j["phantom"]);
    const json& j_R = j["R"];
    PET3D::Matrix<F> R;
    int i = 0;
    for (const auto& el : j_R) {
      R(i++) = el;
    }
    return new RotatedRegion(phantom, R);

  } else if (type == "translated") {
    auto phantom = create_phantom_region_from_json<RNG, F>(j["phantom"]);
    const json& displacement_val = j["displacement"];
    PET3D::Vector<F> displacement;
    displacement.x = displacement_val[0];
    displacement.y = displacement_val[1];
    displacement.z = displacement_val[2];
    return new TranslatedRegion(phantom, displacement);

  } else {
    throw("unknown region type: " + type);
    return nullptr;
  }
}

}  // PET3D

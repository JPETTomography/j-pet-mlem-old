#pragma once

#include "rapidjson/document.h"

#include "point.h"
#include "vector"
#include "event_generator.h"
#include "phantom.h"

namespace PET3D {
using Value = rapidjson::Value;

template <typename AngularDistribution>
AngularDistribution create_angular_distribution_from_json(const Value& obj);

template <>
SingleDirectionDistribution<float> create_angular_distribution_from_json<
    SingleDirectionDistribution<float>>(const Value& obj) {
  using Vector = Vector<float>;

  const Value& type_val = obj["type"];
  if (type_val != "single-direction") {
    std::cerr << "Type mismatch in SingleDirectionDistribution<float> factory: "
              << type_val.GetString() << "\n";
  }

  const rapidjson::Value& direction_val = obj["direction"];
  Vector direction(direction_val[0].GetDouble(),
                   direction_val[1].GetDouble(),
                   direction_val[2].GetDouble());
  return SingleDirectionDistribution<float>(direction);
};

template <>
SphericalDistribution<float> create_angular_distribution_from_json<
    SphericalDistribution<float>>(const Value& obj) {

  const Value& type_val = obj["type"];
  if (type_val != "spherical") {
    std::cerr << "Type mismatch in SphericalDistribution<float> factory: "
              << type_val.GetString() << "\n";
  }
  float theta_min = -M_PI / 2;
  float theta_max = M_PI / 2;
  if (obj.HasMember("theta-min")) {
    theta_min = obj["theta-min"].GetDouble();
  }
  if (obj.HasMember("theta-max")) {
    theta_max = obj["theta-max"].GetDouble();
  }
  return SphericalDistribution<float>(theta_min, theta_max);
}

template <typename FType, typename RNG>
PhantomRegion<FType, RNG>* create_phantom_region_from_json(const Value& obj) {
  if (!obj.HasMember("type")) {
    std::cerr << "phantom region does not have type member\n";
    return nullptr;
  }
  const Value& type_val = obj["type"];
  if (type_val == "cylinder") {
    FType radius = obj["radius"].GetDouble();
    FType height = obj["height"].GetDouble();
    FType intensity = obj["intensity"].GetDouble();

    const Value& angular_val = obj["angular"];
    if (angular_val["type"] == "spherical") {
      SphericalDistribution<FType> angular =
          create_angular_distribution_from_json<SphericalDistribution<FType>>(
              angular_val);
      return new CylinderRegion<FType, RNG>(radius, height, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type_val == "ellipsoid") {
    FType rx = obj["rx"].GetDouble();
    FType ry = obj["ry"].GetDouble();
    FType rz = obj["rz"].GetDouble();
    FType intensity = obj["intensity"].GetDouble();

    const Value& angular_val = obj["angular"];
    if (angular_val["type"] == "spherical") {
      SphericalDistribution<FType> angular =
          create_angular_distribution_from_json<SphericalDistribution<FType>>(
              angular_val);
      return new EllipsoidRegion<FType, RNG>(rx, ry, rz, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  }

  else if (type_val == "point") {
    std::cerr << "not implemented yet\n";
    return nullptr;

  } else if (type_val == "rotated") {
    auto phantom = create_phantom_region_from_json<FType, RNG>(obj["phantom"]);
    const Value& R_val = obj["R"];
    PET3D::Matrix<FType> R;
    int i = 0;
    for (auto it = R_val.Begin(); it != R_val.End(); it++, i++) {
      R(i) = it->GetDouble();
    }

    return new RotatedPhantomRegion<FType, RNG>(phantom, R);

  } else if (type_val == "translated") {

    auto phantom = create_phantom_region_from_json<FType, RNG>(obj["phantom"]);
    const Value& displacement_val = obj["displacement"];
    PET3D::Vector<FType> displacement;
    displacement.x = displacement_val[0].GetDouble();
    displacement.y = displacement_val[1].GetDouble();
    displacement.z = displacement_val[2].GetDouble();

    return new TranslatedPhantomRegion<FType, RNG>(phantom, displacement);

  } else {
    std::cerr << "unknown region type : " << type_val.GetString() << "\n";
    return nullptr;
  }
}
}

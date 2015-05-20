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
PhantomRegion<FType, RNG>* create_phantom_region_from_json(const Value& obj) {}
}

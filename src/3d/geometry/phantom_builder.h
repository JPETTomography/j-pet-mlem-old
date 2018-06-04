#pragma once

#include "point.h"
#include "vector"
#include "event_generator.h"
#include "phantom.h"

#include "util/json.h"

namespace PET3D {

template <class RNG, typename FType>
typename Phantom<RNG, FType>::Region* create_phantom_region_from_json(
    const json& j) {
  using F = FType;
  using Vector = PET3D::Vector<F>;
  using Point = PET3D::Point<F>;
  using Matrix = PET3D::Matrix<F>;
  using SphericalDistribution = Distribution::SphericalDistribution<F>;
  using SingleDirectionDistribution =
      Distribution::SingleDirectionDistribution<F>;
  using Phantom = PET3D::Phantom<RNG, F>;
  using SphericalDistributionPointRegion =
      typename Phantom::template PointRegion<SphericalDistribution>;
  using SingleDirectionDistributionPointRegion =
      typename Phantom::template PointRegion<SingleDirectionDistribution>;
  using CylinderRegion =
      typename Phantom::template CylinderRegion<SphericalDistribution>;
  using EllipsoidRegion =
      typename Phantom::template EllipsoidRegion<SphericalDistribution>;
  using RectangularRegion =
      typename Phantom::template RectangularRegion<SphericalDistribution>;
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
      SphericalDistribution angular(j_angular);
      return new CylinderRegion(radius, height, intensity, angular);
    }

    throw("unsuported cylinder region angular distribution: " + angular_type);

  } else if (type == "ellipsoid") {
    F rx = j["rx"];
    F ry = j["ry"];
    F rz = j["rz"];

    F intensity = j.count("intensity") ? j["intensity"].get<F>() : F(1);

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      SphericalDistribution angular(j_angular);
      return new EllipsoidRegion(rx, ry, rz, intensity, angular);
    }

    throw("unsuported ellipsoid region angular distribution: " + angular_type);

  } else if (type == "point") {
    const json& j_origin = j["origin"];
    Point origin(j_origin[0], j_origin[1], j_origin[2]);
    F intensity = j.count("intensity") ? j["intensity"].get<F>() : F(1);
    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      SphericalDistribution angular(j_angular);
      return new SphericalDistributionPointRegion(intensity, angular, origin);
    } else if (angular_type == "single-direction") {
      SingleDirectionDistribution angular(j_angular);
      return new SingleDirectionDistributionPointRegion(
          intensity, angular, origin);
    }

    throw("unsupported point region angular distribution: " + angular_type);

  } else if (type == "rectangular") {
    F ax = j["ax"];
    F ay = j["ay"];
    F az = j["az"];
    F intensity = j.count("intensity") ? j["intensity"].get<F>() : F(1);
    return new RectangularRegion(ax, ay, az, intensity);

  } else if (type == "rotated") {
    auto phantom = create_phantom_region_from_json<RNG, F>(j["phantom"]);
    const json& j_R = j["R"];
    const json& j_axis = j["axis"];
    const json& j_angle = j["angle"];
    if (j_R.is_array()) {
      Matrix R;
      int i = 0;
      for (const auto& el : j_R) {
        R[i++] = el;
      }
      return new RotatedRegion(phantom, R);
    } else if (j_axis.is_array() && j_angle.is_number()) {
      auto R = Matrix::rotate(Vector(j_axis[0], j_axis[1], j_axis[2]),
                              (double)j_angle * M_PI / 180);
      return new RotatedRegion(phantom, R);
    }

    throw("rotated phantom must contain axis & angle pair or R matrix");

  } else if (type == "translated") {
    auto phantom = create_phantom_region_from_json<RNG, F>(j["phantom"]);
    const json& displacement_val = j["displacement"];
    Vector displacement;
    displacement.x = displacement_val[0];
    displacement.y = displacement_val[1];
    displacement.z = displacement_val[2];
    return new TranslatedRegion(phantom, displacement);
  }

  throw("unknown region type: " + type);
}

}  // PET3D

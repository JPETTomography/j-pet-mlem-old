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
  using Matrix = PET3D::Matrix<F>;
  using AngularDistribution = Distribution::SphericalDistribution<F>;
  using Phantom = PET3D::Phantom<RNG, F>;
  using CylinderRegion =
      typename Phantom::template CylinderRegion<AngularDistribution>;
  using EllipsoidRegion =
      typename Phantom::template EllipsoidRegion<AngularDistribution>;
  using RectangularRegion =
      typename Phantom::template RectangularRegion<AngularDistribution>;
  using RotatedRegion = typename Phantom::RotatedRegion;
  using TranslatedRegion = typename Phantom::TranslatedRegion;

  if (!j.count("type")) {
    throw("phantom region does not have type member");
  }

  std::string type = j["type"];
  std::cerr << "type : " << type << "\n";

  if (type == "cylinder") {
    F radius = j["radius"];
    F height = j["height"];
    F intensity = j["intensity"];

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      AngularDistribution angular(j_angular);
      return new CylinderRegion(radius, height, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "ellipsoid") {
    std::cerr << "ellipsoid ";
    F rx = j["rx"];
    F ry = j["ry"];
    F rz = j["rz"];
    std::cerr << rx << " " << ry << " " << rz << "\n";

    F intensity = j.count("intensity") ? j["intensity"].get<F>() : F(1);

    const json& j_angular = j["angular"];
    std::string angular_type = j_angular["type"];
    if (angular_type == "spherical") {
      AngularDistribution angular(j_angular);
      return new EllipsoidRegion(rx, ry, rz, intensity, angular);
    } else {
      std::cerr << "unsuported angular distribution\n";
      return nullptr;
    }

  } else if (type == "point") {
    throw("not implemented yet");

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
    } else {
      throw("rotated phantom must contain axis & angle pair or R matrix");
    }

  } else if (type == "translated") {
    auto phantom = create_phantom_region_from_json<RNG, F>(j["phantom"]);
    const json& displacement_val = j["displacement"];
    Vector displacement;
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

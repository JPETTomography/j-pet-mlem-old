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
SingleDirectionDistribution<float>
create_angular_distribution_from_json<SingleDirectionDistribution<float>>(const Value& obj) {
  using Vector = Vector<float>;

  const Value& type_val = obj["type"];
  if(type_val != "single-direction") {
      std::cerr<<"Type mismatch in SingleDirectionDistribution<float> factory: "<<type_val.GetString()<<"\n";
  }

  const rapidjson::Value& direction_val = obj["direction"];
  Vector direction(direction_val[0].GetDouble(), direction_val[1].GetDouble(), direction_val[2].GetDouble() );
  return SingleDirectionDistribution<float>(direction);
};



template<typename FType, typename RNG>
PhantomRegion<FType, RNG>* create_phantom_region_from_json(const Value& obj) {

}


}

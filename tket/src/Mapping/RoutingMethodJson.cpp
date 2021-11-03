#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Utils/Json.hpp"

namespace tket{

void to_json(nlohmann::json& j, const RoutingMethod& rm) { j = rm.serialize(); }

void from_json(nlohmann::json& j, RoutingMethod& rm) {
  std::string name = j.at("name").get<std::string>();
  if (name == "LexiRouteRoutingMethod") {
      rm = LexiRouteRoutingMethod::deserialize(j);
  } else {
    throw JsonError(
        "Deserialization not yet implemented for generic RoutingMethod "
        "objects.");
  }
}

JSON_DECL(RoutingMethod);


void to_json(nlohmann::json& j, const
std::vector<RoutingMethodWrapper>& r){
  std::vector<nlohmann::json> serialized_methods;  
  for(const RoutingMethodWrapper& routing_method : r){
    serialized_methods.push_back(routing_method.get().serialize());
  }
  j["methods"] = serialized_methods;
}

void from_json(nlohmann::json& j,
std::vector<RoutingMethodWrapper>& r){
  std::vector<nlohmann::json> serialized_methods = j.at("methods").get<std::vector<nlohmann::json>>();
  for(const nlohmann::json& nj : serialized_methods){
      std::string name = nj.at("name").get<std::string>();
      if(name == "LexiRouteRoutingMethod"){
          LexiRouteRoutingMethod lrrm = LexiRouteRoutingMethod::deserialize(nj);
          r.push_back(lrrm);
      }
      else{
          throw JsonError("Deserialization not implemented for generic RoutingMethod objects in vector.");
      }
  }
}
}

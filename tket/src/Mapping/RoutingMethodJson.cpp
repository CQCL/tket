#include "Mapping/RoutingMethodJson.hpp"

namespace tket {

void to_json(nlohmann::json& j, const RoutingMethod& rm) { j = rm.serialize(); }

void from_json(const nlohmann::json& j, RoutingMethod& rm) {
  std::string name = j.at("name").get<std::string>();
  if (name == "LexiRouteRoutingMethod") {
    rm = LexiRouteRoutingMethod::deserialize(j);
  } else {
    throw JsonError(
        "Deserialization not yet implemented for generic RoutingMethod "
        "objects.");
  }
}

}  // namespace tket

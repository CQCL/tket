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

void to_json(nlohmann::json& j, const std::vector<RoutingMethodPtr>& rmp) {
  for (const auto& r : rmp) {
    j.push_back(*r);
  }
}

void from_json(const nlohmann::json& j, std::vector<RoutingMethodPtr>& rmp) {
  for (const auto& c : j) {
    std::string name = c.at("name").get<std::string>();
    if (name == "LexiRouteRoutingMethod") {
      rmp.push_back(std::make_shared<LexiRouteRoutingMethod>(
          LexiRouteRoutingMethod::deserialize(c)));
    } else {
      rmp.push_back(std::make_shared<RoutingMethod>(c.get<RoutingMethod>()));
    }
  }
}

}  // namespace tket

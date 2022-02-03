#include "Mapping/RoutingMethodJson.hpp"

namespace tket {

void to_json(nlohmann::json& j, const RoutingMethod& rm) { j = rm.serialize(); }

void from_json(const nlohmann::json& /*j*/, RoutingMethod& rm) {
  rm = RoutingMethod();
}

void to_json(nlohmann::json& j, const std::vector<RoutingMethodPtr>& rmp_v) {
  for (const auto& r : rmp_v) {
    j.push_back(*r);
  }
}

void from_json(const nlohmann::json& j, std::vector<RoutingMethodPtr>& rmp_v) {
  for (const auto& c : j) {
    std::string name = c.at("name").get<std::string>();
    if (name == "LexiRouteRoutingMethod") {
      rmp_v.push_back(std::make_shared<LexiRouteRoutingMethod>(
          LexiRouteRoutingMethod::deserialize(c)));
    } else if (name == "RoutingMethod") {
      rmp_v.push_back(std::make_shared<RoutingMethod>());
    } else if (name == "MultiGateReorderRoutingMethod") {
      rmp_v.push_back(std::make_shared<MultiGateReorderRoutingMethod>(
          MultiGateReorderRoutingMethod::deserialize(c)));
    } else {
      std::logic_error(
          "Deserialization for given RoutingMethod not supported.");
    }
  }
}

}  // namespace tket

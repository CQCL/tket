#ifndef _TKET_RoutingMethodJson_H_
#define _TKET_RoutingMethodJson_H_

#include "Mapping/LexiRoute.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Utils/Json.hpp"

namespace tket {

void to_json(nlohmann::json& j, const RoutingMethod& rm);

void from_json(const nlohmann::json& j, RoutingMethod& rm);

JSON_DECL(RoutingMethod);

void to_json(nlohmann::json& j, const std::vector<RoutingMethodPtr>& rmp);

void from_json(const nlohmann::json& j, std::vector<RoutingMethodPtr>& rmp);

JSON_DECL(std::vector<RoutingMethodPtr>);

}  // namespace tket

#endif

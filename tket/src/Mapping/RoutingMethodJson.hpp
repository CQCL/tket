#ifndef _TKET_RoutingMethodJson_H_
#define _TKET_RoutingMethodJson_H_

#include "Mapping/LexiRoute.hpp"
#include "Utils/Json.hpp"

namespace tket{

void to_json(nlohmann::json& j, const RoutingMethod& rm);

void from_json(nlohmann::json& j, RoutingMethod& rm);

JSON_DECL(RoutingMethod);

void to_json(nlohmann::json& j, const
std::vector<RoutingMethodWrapper>& r);

void from_json(nlohmann::json& j,
std::vector<RoutingMethodWrapper>& r);

JSON_DECL(std::vector<RoutingMethodWrapper>);

} // namespace tket

#endif

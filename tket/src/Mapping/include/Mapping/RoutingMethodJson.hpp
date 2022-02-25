// Copyright 2019-2022 Cambridge Quantum Computing
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include "Mapping/AASLabelling.hpp"
#include "Mapping/AASRoute.hpp"
#include "Mapping/BoxDecomposition.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRouteRoutingMethod.hpp"
#include "Mapping/MultiGateReorder.hpp"
#include "Mapping/RoutingMethod.hpp"
#include "Utils/Json.hpp"

namespace tket {

void to_json(nlohmann::json& j, const RoutingMethod& rm);

void from_json(const nlohmann::json& /*j*/, RoutingMethod& rm);

JSON_DECL(RoutingMethod);

void to_json(nlohmann::json& j, const std::vector<RoutingMethodPtr>& rmp_v);

void from_json(const nlohmann::json& j, std::vector<RoutingMethodPtr>& rmp_v);

JSON_DECL(std::vector<RoutingMethodPtr>);

}  // namespace tket

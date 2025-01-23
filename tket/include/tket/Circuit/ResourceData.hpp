// Copyright Quantinuum
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

#include <concepts>
#include <map>

#include "tket/OpType/OpType.hpp"

namespace tket {

template <typename T>
concept arithmetic = std::integral<T> or std::floating_point<T>;

template <arithmetic T>
struct ResourceBounds {
  T min;
  T max;
  bool operator==(const ResourceBounds<T>& other) const {
    return min == other.min && max == other.max;
  }
  ResourceBounds() : min(0), max(0) {}
  ResourceBounds(T val) : min(val), max(val) {}
  ResourceBounds(T minval, T maxval) : min(minval), max(maxval) {}
};

struct ResourceData {
  std::map<OpType, ResourceBounds<unsigned>> OpTypeCount;
  ResourceBounds<unsigned> GateDepth;
  std::map<OpType, ResourceBounds<unsigned>> OpTypeDepth;
  ResourceBounds<unsigned> TwoQubitGateDepth;
  bool operator==(const ResourceData& other) const;
};

template <arithmetic T>
void to_json(nlohmann::json& j, const ResourceBounds<T>& bounds) {
  j["min"] = bounds.min;
  j["max"] = bounds.max;
}
template <arithmetic T>
void from_json(const nlohmann::json& j, ResourceBounds<T>& bounds) {
  bounds.min = j.at("min").get<T>();
  bounds.max = j.at("max").get<T>();
}

JSON_DECL(ResourceData)

}  // namespace tket

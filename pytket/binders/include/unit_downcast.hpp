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
#include <nanobind/nanobind.h>

#include <typeinfo>

#include "tket/Utils/UnitID.hpp"

namespace nanobind::detail {

/** Enable automatic downcasting of UnitIDs as required for some Circuit methods
 */
template <>
struct type_hook<tket::UnitID> {
  static const std::type_info *get(tket::UnitID *src) {
    if (src) {
      if (src->type() == tket::UnitType::Qubit) {
        // Node has no additional info but is more specific
        // If Qubit is needed, then subtyping is sufficient
        return &typeid(tket::Node);
      } else if (src->type() == tket::UnitType::WasmState) {
        return &typeid(tket::WasmState);
      } else {
        return &typeid(tket::Bit);
      }
    }
    return &typeid(tket::UnitID);
  }
};
}  // namespace nanobind::detail

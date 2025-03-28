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

#include "tket/Utils/UnitID.hpp"

namespace nanobind {

// ((( Copied from pybind11 (type_caster_base.h)
template <typename itype, typename SFINAE = void>
struct polymorphic_type_hook_base {
    static const void *get(const itype *src, const std::type_info *&) { return src; }
};
template <typename itype>
struct polymorphic_type_hook_base<itype, detail::enable_if_t<std::is_polymorphic<itype>::value>> {
    static const void *get(const itype *src, const std::type_info *&type) {
        type = src ? &typeid(*src) : nullptr;
        return dynamic_cast<const void *>(src);
    }
};
template <typename itype, typename SFINAE = void>
struct polymorphic_type_hook : public polymorphic_type_hook_base<itype> {};
// )))

/** Enable automatic downcasting of UnitIDs as required for some Circuit methods
 */
template <>
struct polymorphic_type_hook<tket::UnitID> {
  static const void* get(const tket::UnitID* src, const std::type_info*& type) {
    if (src) {
      if (src->type() == tket::UnitType::Qubit) {
        // Node has no additional info but is more specific
        // If Qubit is needed, then subtyping is sufficient
        type = &typeid(tket::Node);
      } else if (src->type() == tket::UnitType::WasmState) {
        type = &typeid(tket::WasmState);
      } else {
        type = &typeid(tket::Bit);
      }
    } else
      type = nullptr;
    return src;
  }
};
}  // namespace nanobind

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

#include "UnitRegister.hpp"
#include "tket/Circuit/Circuit.hpp"

namespace tket {

template <typename T>
std::vector<T> get_unit_registers(Circuit &circ) {
  static_assert(
      std::is_same<T, QubitRegister>::value ||
          std::is_same<T, BitRegister>::value,
      "T must be either QubitRegister or BitRegister");
  using T2 = typename std::conditional<
      std::is_same<T, QubitRegister>::value, Qubit, Bit>::type;
  std::vector<T2> unitids;
  if constexpr (std::is_same<T, QubitRegister>::value) {
    unitids = circ.all_qubits();
  } else if constexpr (std::is_same<T, BitRegister>::value) {
    unitids = circ.all_bits();
  }
  // map from register name to unsigned indices
  std::map<std::string, std::set<unsigned>> unit_map;
  std::vector<T> regs;
  for (const T2 &unitid : unitids) {
    // UnitRegisters only describe registers with 1-d indices
    if (unitid.reg_dim() != 1) continue;
    auto it = unit_map.find(unitid.reg_name());
    if (it == unit_map.end()) {
      unit_map.insert({unitid.reg_name(), {unitid.index()[0]}});
    } else {
      it->second.insert(unitid.index()[0]);
    }
  }
  regs.reserve(unit_map.size());
  for (auto const &it : unit_map) {
    // only return registers that are indexed consecutively from zero
    if (*it.second.rbegin() == it.second.size() - 1) {
      regs.emplace_back(it.first, it.second.size());
    }
  }
  return regs;
}

}  // namespace tket

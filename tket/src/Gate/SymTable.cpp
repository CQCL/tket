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

#include "SymTable.hpp"

namespace tket {

std::unordered_set<std::string>& SymTable::get_registered_symbols() {
  static std::unordered_set<std::string> symbols;
  return symbols;
}

Sym SymTable::fresh_symbol(const std::string& preferred) {
  std::string new_symbol = preferred;
  unsigned suffix = 0;
  while (get_registered_symbols().find(new_symbol) !=
         get_registered_symbols().cend()) {
    suffix++;
    new_symbol = preferred + "_" + std::to_string(suffix);
  }
  register_symbol(new_symbol);
  return SymEngine::symbol(new_symbol);
}

void SymTable::register_symbol(const std::string& symbol) {
  get_registered_symbols().insert(symbol);
}

void SymTable::register_symbols(const SymSet& ss) {
  for (const auto& s : ss) {
    get_registered_symbols().insert(s->get_name());
  }
}

}  // namespace tket

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

#include "DiagUtils.hpp"

namespace tket {

void insert_into_gadget_map(
    QubitOperator &gadget_map, const PauliGadgetProperties &pgp) {
  QubitOperator::iterator iter = gadget_map.find(pgp.tensor_);
  if (iter == gadget_map.end())
    gadget_map[pgp.tensor_] = pgp.angle_;
  else {
    QubitPauliTensor string_to_insert = pgp.tensor_ * iter->first;
    Expr ang_to_insert = pgp.angle_ * iter->second;
    gadget_map.erase(iter);
    gadget_map[string_to_insert] = ang_to_insert;
  }
}

}  // namespace tket

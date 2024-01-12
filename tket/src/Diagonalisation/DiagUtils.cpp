// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tket/Diagonalisation/DiagUtils.hpp"

namespace tket {

void insert_into_gadget_map(
    QubitOperator &gadget_map, const PauliGadgetProperties &pgp) {
  SpSymPauliTensor gadget(pgp.tensor_);
  gadget.coeff *= pgp.angle_;
  SpPauliString ps(gadget);
  QubitOperator::iterator iter = gadget_map.find(ps);
  if (iter == gadget_map.end())
    gadget_map[ps] = gadget.coeff;
  else {
    iter->second += gadget.coeff;
  }
}

}  // namespace tket

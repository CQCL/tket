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

#include "tket/Converters/Converters.hpp"

namespace tket {

APState circuit_to_apstate(const Circuit& circ) {
  APState aps(circ.n_qubits());
  std::map<UnitID, unsigned> qb_ordering;
  for (const Qubit& q : circ.all_qubits())
    qb_ordering.insert({q, qb_ordering.size()});
  for (const Command& com : circ) {
    auto args = com.get_args();
    std::vector<unsigned> qbs;
    for (const UnitID& q : args) qbs.push_back(qb_ordering.at(q));
    aps.apply_gate(com.get_op_ptr()->get_type(), qbs);
  }
  return aps;
}

}  // namespace tket

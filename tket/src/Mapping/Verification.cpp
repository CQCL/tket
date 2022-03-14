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

#include "Verification.hpp"

namespace tket {
bool respects_connectivity_constraints(
    const Circuit &circ, const Architecture &arch, bool directed,
    bool bridge_allowed) {
  unit_set_t qb_lookup;
  for (const Qubit &qb : circ.all_qubits()) {
    if (!arch.node_exists(Node(qb))) return false;
    qb_lookup.insert(qb);
  }
  for (Circuit::CommandIterator it = circ.begin(); it != circ.end(); ++it) {
    const Command &com = *it;
    unit_vector_t qbs;
    for (const UnitID &arg : com.get_args()) {
      if (qb_lookup.find(arg) != qb_lookup.end()) qbs.push_back(arg);
    }
    Op_ptr op = com.get_op_ptr();
    if (op->get_type() == OpType::Barrier) continue;
    if (op->get_type() == OpType::Conditional) {
      std::shared_ptr<const Conditional> cond_ptr =
          std::dynamic_pointer_cast<const Conditional>(op);
      op = cond_ptr->get_op();
    }
    if (op->get_type() == OpType::CircBox) {
      std::shared_ptr<const Box> box_ptr =
          std::dynamic_pointer_cast<const Box>(op);
      Circuit box_circ = *box_ptr->to_circuit().get();
      qubit_vector_t all_units = box_circ.all_qubits();
      if (all_units.size() != qbs.size()) return false;
      unit_map_t rename_map;
      for (unsigned i = 0; i < all_units.size(); i++)
        rename_map.insert({all_units[i], qbs[i]});
      box_circ.rename_units(rename_map);
      if (!respects_connectivity_constraints(
              box_circ, arch, directed, bridge_allowed))
        return false;
      continue;
    }
    unsigned n_qbs = qbs.size();
    switch (n_qbs) {
      case 1:
        break;
      case 2: {
        if (arch.get_distance(Node(qbs[0]), Node(qbs[1])) != 1) {
          return false;
        }
        if (directed) {
          OpType ot = op->get_type();
          if ((ot == OpType::CX || ot == OpType::ECR) &&
              !arch.edge_exists(Node(qbs[0]), Node(qbs[1])))
            return false;
        }
        break;
      }
      case 3: {
        if (bridge_allowed) {
          if (directed)
            throw std::logic_error(
                "BRIDGE ops are disallowed on a directed "
                "architecture. They must be decomposed.");
          if (op->get_type() == OpType::BRIDGE) {
            if (arch.get_distance(Node(qbs[0]), Node(qbs[1])) != 1 ||
                arch.get_distance(Node(qbs[1]), Node(qbs[2])) != 1) {
              return false;
            }
          } else
            return false;
        } else
          return false;
        break;
      }
      default:
        return false;
    }
  }
  return true;
}

}  // namespace tket

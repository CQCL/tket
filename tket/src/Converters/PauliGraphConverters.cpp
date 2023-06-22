// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/Converters/PauliGadget.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Diagonalisation/DiagUtils.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"
#include "tket/Gate/Gate.hpp"

namespace tket {

PauliGraph circuit_to_pauli_graph(const Circuit &circ) {
  PauliGraph pg(circ.all_qubits(), circ.all_bits());
  for (const Command &com : circ) {
    const Op &op = *com.get_op_ptr();
    unit_vector_t args = com.get_args();
    OpDesc od = op.get_desc();
    if (od.is_gate()) {
      pg.apply_gate_at_end(static_cast<const Gate &>(op), args);
    } else if (od.type() == OpType::PauliExpBox) {
      const PauliExpBox &peb = static_cast<const PauliExpBox &>(op);
      std::vector<Pauli> paulis = peb.get_paulis();
      Expr phase = peb.get_phase();
      if (args.size() != paulis.size())
        throw std::logic_error("Incorrect Pauli tensor size for qubit count");
      QubitPauliMap qpm;
      for (unsigned i = 0; i != args.size(); ++i)
        qpm.insert({Qubit(args[i]), paulis[i]});
      QubitPauliTensor qpt = pg.cliff_.get_row_product(QubitPauliTensor(qpm));
      pg.apply_pauli_gadget_at_end(qpt, phase);
    } else
      throw BadOpType(
          "Can only make a PauliGraph from a circuit of basic gates "
          "and Paulis",
          od.type());
  }
  return pg;
}

Circuit pauli_graph_to_pauli_exp_box_circuit_individually(
    const PauliGraph &pg, CXConfigType cx_config) {
  Circuit circ;
  for (const Qubit &qb : pg.cliff_.get_qubits()) {
    circ.add_qubit(qb);
  }
  for (const Bit &b : pg.bits_) {
    circ.add_bit(b);
  }
  for (const PauliVert &vert : pg.vertices_in_order()) {
    const QubitPauliTensor &pauli = pg.graph_[vert].tensor_;
    const Expr &angle = pg.graph_[vert].angle_;
    append_single_pauli_gadget_as_pauli_exp_box(circ, pauli, angle, cx_config);
  }
  Circuit cliff_circuit = unitary_rev_tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it = pg.measures_.begin(); it != pg.measures_.end(); ++it) {
    circ.add_measure(it->left, it->right);
  }
  return circ;
}

Circuit pauli_graph_to_pauli_exp_box_circuit_pairwise(
    const PauliGraph &pg, CXConfigType cx_config) {
  Circuit circ;
  for (const Qubit &qb : pg.cliff_.get_qubits()) {
    circ.add_qubit(qb);
  }
  for (const Bit &b : pg.bits_) {
    circ.add_bit(b);
  }
  std::vector<PauliVert> vertices = pg.vertices_in_order();
  auto it = vertices.begin();
  while (it != vertices.end()) {
    PauliVert vert0 = *it;
    const QubitPauliTensor &pauli0 = pg.graph_[vert0].tensor_;
    const Expr &angle0 = pg.graph_[vert0].angle_;
    ++it;
    if (it == vertices.end()) {
      // append_single_pauli_gadget(circ, pauli0, angle0, cx_config);
      append_single_pauli_gadget_as_pauli_exp_box(
          circ, pauli0, angle0, cx_config);
    } else {
      PauliVert vert1 = *it;
      const QubitPauliTensor &pauli1 = pg.graph_[vert1].tensor_;
      const Expr &angle1 = pg.graph_[vert1].angle_;
      ++it;
      append_pauli_gadget_pair_as_box(
          circ, pauli0, angle0, pauli1, angle1, cx_config);
      // append_pauli_gadget_pair(circ, pauli0, angle0, pauli1, angle1,
      // cx_config);
    }
  }
  Circuit cliff_circuit = unitary_rev_tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it = pg.measures_.begin(); it != pg.measures_.end(); ++it) {
    circ.add_measure(it->left, it->right);
  }
  return circ;
}

/**
 * @brief Greedily group gadgets into commuting sets
 *
 * @param pg
 * @param graph so we don't need to make this function a friend
 * @return std::vector<QubitOperator>
 */
static std::vector<QubitOperator> group_commuting_gagdets(
    const PauliGraph &pg, const PauliDAG &graph) {
  std::vector<QubitOperator> commuting_gagdets;
  std::vector<PauliVert> vertices = pg.vertices_in_order();
  auto it = vertices.begin();
  while (it != vertices.end()) {
    const PauliGadgetProperties &pgp = graph[*it];
    QubitOperator gadget_map;
    gadget_map[pgp.tensor_] = pgp.angle_;
    ++it;
    while (it != vertices.end()) {
      const PauliGadgetProperties &pauli_gadget = graph[*it];
      QubitOperator::iterator pgs_iter = gadget_map.find(pauli_gadget.tensor_);
      if (pgs_iter != gadget_map.end()) {
        insert_into_gadget_map(gadget_map, pauli_gadget);
      } else {
        bool commutes_with_all = true;
        for (const std::pair<const QubitPauliTensor, Expr> &pv : gadget_map) {
          if (!pauli_gadget.tensor_.commutes_with(pv.first)) {
            commutes_with_all = false;
            break;
          }
        }
        if (!commutes_with_all) break;
        insert_into_gadget_map(gadget_map, pauli_gadget);
      }
      ++it;
    }
    commuting_gagdets.push_back(gadget_map);
  }
  return commuting_gagdets;
}

///* Currently follows a greedy set-building method */
Circuit pauli_graph_to_pauli_exp_box_circuit_sets(
    const PauliGraph &pg, CXConfigType cx_config) {
  Circuit circ;
  const std::set<Qubit> qbs = pg.cliff_.get_qubits();
  Circuit spare_circ;
  for (const Qubit &qb : qbs) {
    circ.add_qubit(qb);
    spare_circ.add_qubit(qb);
  }
  for (const Bit &b : pg.bits_) {
    circ.add_bit(b);
  }
  std::vector<QubitOperator> commuting_gagdets =
      group_commuting_gagdets(pg, pg.graph_);
  for (const QubitOperator &gadget_map : commuting_gagdets) {
    if (gadget_map.size() == 1) {
      const std::pair<const QubitPauliTensor, Expr> &pgp0 = *gadget_map.begin();
      append_single_pauli_gadget_as_pauli_exp_box(
          circ, pgp0.first, pgp0.second, cx_config);
    } else if (gadget_map.size() == 2) {
      const std::pair<const QubitPauliTensor, Expr> &pgp0 = *gadget_map.begin();
      const std::pair<const QubitPauliTensor, Expr> &pgp1 =
          *(++gadget_map.begin());
      append_pauli_gadget_pair_as_box(
          circ, pgp0.first, pgp0.second, pgp1.first, pgp1.second, cx_config);
    } else {
      std::list<std::pair<QubitPauliTensor, Expr>> gadgets;
      for (const std::pair<const QubitPauliTensor, Expr> &qps_pair :
           gadget_map) {
        gadgets.push_back(qps_pair);
      }
      append_commuting_pauli_gadget_set_as_box(circ, gadgets, cx_config);
    }
  }
  Circuit cliff_circuit = unitary_rev_tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it1 = pg.measures_.begin(); it1 != pg.measures_.end(); ++it1) {
    circ.add_measure(it1->left, it1->right);
  }
  return circ;
}

}  // namespace tket

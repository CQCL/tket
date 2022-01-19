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

#include "Circuit/Boxes.hpp"
#include "Converters.hpp"
#include "Converters/PhasePoly.hpp"
#include "Diagonalisation/Diagonalisation.hpp"
#include "Gate/Gate.hpp"
#include "PauliGadget.hpp"

namespace tket {

PauliGraph circuit_to_pauli_graph(const Circuit &circ);
Circuit pauli_graph_to_circuit(const PauliGraph &pg);

PauliGraph circuit_to_pauli_graph(const Circuit &circ) {
  for (const std::pair<const Qubit, Qubit> &pair :
       circ.implicit_qubit_permutation()) {
    if (pair.first != pair.second) {
      throw NotImplemented(
          "Cannot build a PauliGraph from circuits with implicit "
          "permutations");
    }
  }
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
        throw NotValid("Incorrect Pauli tensor size for qubit count");
      QubitPauliTensor qpt;
      for (unsigned i = 0; i != args.size(); ++i) {
        switch (paulis[i]) {
          case Pauli::I:
            break;
          case Pauli::X:
            qpt = qpt * pg.cliff_.get_xpauli(Qubit(args[i]));
            break;
          case Pauli::Y:
            qpt = qpt * pg.cliff_.get_xpauli(Qubit(args[i]));
            qpt = qpt * pg.cliff_.get_zpauli(Qubit(args[i]));
            qpt = i_ * qpt;
            break;
          case Pauli::Z:
            qpt = qpt * pg.cliff_.get_zpauli(Qubit(args[i]));
            break;
        }
      }
      pg.apply_pauli_gadget_at_end(qpt, phase);
    } else
      throw NotImplemented(
          "Can only make a PauliGraph from a circuit of basic gates "
          "and Paulis");
  }
  return pg;
}

Circuit pauli_graph_to_circuit_individually(
    const PauliGraph &pg, CXConfigType cx_config) {
  Circuit circ;
  for (const Qubit &qb : pg.cliff_.get_qubits()) {
    circ.add_qubit(qb);
  }
  for (const Bit &b : pg.bits_) {
    circ.add_bit(b);
  }
  for (PauliGraph::TopSortIterator it = pg.begin(); it != pg.end(); ++it) {
    PauliVert vert = *it;
    const QubitPauliTensor &pauli = pg.graph_[vert].tensor_;
    const Expr &angle = pg.graph_[vert].angle_;
    append_single_pauli_gadget(circ, pauli, angle, cx_config);
  }
  Circuit cliff_circuit = tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it = pg.measures_.begin(); it != pg.measures_.end(); ++it) {
    circ.add_measure(it->left, it->right);
  }
  return circ;
}

Circuit pauli_graph_to_circuit_pairwise(
    const PauliGraph &pg, CXConfigType cx_config) {
  Circuit circ;
  for (const Qubit &qb : pg.cliff_.get_qubits()) {
    circ.add_qubit(qb);
  }
  for (const Bit &b : pg.bits_) {
    circ.add_bit(b);
  }
  PauliGraph::TopSortIterator it = pg.begin();
  while (it != pg.end()) {
    PauliVert vert0 = *it;
    const QubitPauliTensor &pauli0 = pg.graph_[vert0].tensor_;
    const Expr &angle0 = pg.graph_[vert0].angle_;
    ++it;
    if (it == pg.end()) {
      append_single_pauli_gadget(circ, pauli0, angle0, cx_config);
    } else {
      PauliVert vert1 = *it;
      const QubitPauliTensor &pauli1 = pg.graph_[vert1].tensor_;
      const Expr &angle1 = pg.graph_[vert1].angle_;
      ++it;
      append_pauli_gadget_pair(circ, pauli0, angle0, pauli1, angle1, cx_config);
    }
  }
  Circuit cliff_circuit = tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it = pg.measures_.begin(); it != pg.measures_.end(); ++it) {
    circ.add_measure(it->left, it->right);
  }
  return circ;
}

/* Currently follows a greedy set-building method */
Circuit pauli_graph_to_circuit_sets(
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
  PauliGraph::TopSortIterator it = pg.begin();
  while (it != pg.end()) {
    const PauliGadgetProperties &pgp = pg.graph_[*it];
    QubitOperator gadget_map;
    gadget_map[pgp.tensor_] = pgp.angle_;
    ++it;
    while (it != pg.end()) {
      const PauliGadgetProperties &pauli_gadget = pg.graph_[*it];
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
    if (gadget_map.size() == 1) {
      const std::pair<const QubitPauliTensor, Expr> &pgp0 = *gadget_map.begin();
      append_single_pauli_gadget(circ, pgp0.first, pgp0.second, cx_config);
    } else if (gadget_map.size() == 2) {
      const std::pair<const QubitPauliTensor, Expr> &pgp0 = *gadget_map.begin();
      const std::pair<const QubitPauliTensor, Expr> &pgp1 =
          *(++gadget_map.begin());
      append_pauli_gadget_pair(
          circ, pgp0.first, pgp0.second, pgp1.first, pgp1.second, cx_config);
    } else {
      std::list<std::pair<QubitPauliTensor, Expr>> gadgets;
      for (const std::pair<const QubitPauliTensor, Expr> &qps_pair :
           gadget_map) {
        gadgets.push_back(qps_pair);
      }
      Circuit cliff_circ = mutual_diagonalise(gadgets, qbs, cx_config);
      circ.append(cliff_circ);
      Circuit phase_poly_circ(spare_circ);
      for (const std::pair<QubitPauliTensor, Expr> &pgp : gadgets) {
        append_single_pauli_gadget(phase_poly_circ, pgp.first, pgp.second);
      }
      PhasePolyBox ppbox(phase_poly_circ);
      Circuit after_synth_circ = *ppbox.to_circuit();
      circ.append(after_synth_circ);
      circ.append(cliff_circ.dagger());
    }
  }
  Circuit cliff_circuit = tableau_to_circuit(pg.cliff_);
  circ.append(cliff_circuit);
  for (auto it = pg.measures_.begin(); it != pg.measures_.end(); ++it) {
    circ.add_measure(it->left, it->right);
  }
  return circ;
}

}  // namespace tket

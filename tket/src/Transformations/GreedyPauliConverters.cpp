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

#include <algorithm>

#include "tket/Circuit/PauliExpBoxes.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/PauliGraph/PauliGraph.hpp"
#include "tket/Transformations/CliffordOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisation.hpp"
#include "tket/Transformations/GreedyPauliOptimisationLookupTables.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

namespace GreedyPauliSimp {

// convert a Pauli exponential to a PauliNode_ptr
static PauliNode_ptr get_node_from_exp(
    const std::vector<Pauli>& paulis, const Expr& theta,
    const qubit_vector_t& args, unsigned n) {
  // pad the Paulis
  std::vector<Pauli> string(n, Pauli::I);
  for (unsigned i = 0; i < args.size(); i++) {
    string[args[i].index()[0]] = paulis[i];
  }
  return std::make_shared<PauliRotation>(string, theta);
}

// convert a Clifford tableau to a vector of PauliNode_ptr
static std::vector<PauliNode_ptr> get_nodes_from_tableau(
    const UnitaryRevTableau& tab, unsigned n_qubits) {
  std::vector<PauliNode_ptr> rows;
  for (unsigned i = 0; i < n_qubits; i++) {
    Qubit q(i);
    SpPauliStabiliser z_stab = tab.get_zrow(q);
    SpPauliStabiliser x_stab = tab.get_xrow(q);
    bool z_sign = cast_coeff<quarter_turns_t, Complex>(z_stab.coeff) == 1.;
    bool x_sign = cast_coeff<quarter_turns_t, Complex>(x_stab.coeff) == 1.;
    TKET_ASSERT(z_stab.string.size() == n_qubits);
    std::vector<Pauli> z_string;
    std::vector<Pauli> x_string;
    for (unsigned j = 0; j < n_qubits; j++) {
      z_string.push_back(z_stab.string.at(Qubit(j)));
      x_string.push_back(x_stab.string.at(Qubit(j)));
    }
    rows.push_back(std::make_shared<PauliPropagation>(
        z_string, x_string, z_sign, x_sign, i));
  }
  return rows;
}

// detect trivial pauli exps, if true then return the global phase
static std::pair<bool, Expr> is_trivial_pauliexp(
    const std::vector<Pauli>& paulis, const Expr& theta) {
  if (static_cast<std::size_t>(std::count(
          paulis.begin(), paulis.end(), Pauli::I)) == paulis.size()) {
    // If all identity term
    return {true, -theta / 2};
  }
  if (equiv_0(theta, 2)) {
    if (equiv_0(theta, 4)) {
      return {true, 0};
    } else {
      return {true, -1};
    }
  }
  return {false, 0};
}

std::tuple<std::vector<PauliNode_ptr>, std::vector<PauliNode_ptr>>
gpg_from_unordered_set(const std::vector<SymPauliTensor>& unordered_set) {
  std::vector<PauliNode_ptr> rotation_set;
  unsigned n_qubits = unordered_set[0].string.size();
  for (auto& pauli : unordered_set) {
    TKET_ASSERT(pauli.string.size() == n_qubits);
    rotation_set.push_back(
        std::make_shared<PauliRotation>(pauli.string, pauli.coeff));
  }
  UnitaryRevTableau tab(n_qubits);
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  return {rotation_set, rows};
}

std::tuple<
    Circuit, std::vector<std::vector<PauliNode_ptr>>,
    std::vector<PauliNode_ptr>, Circuit, unit_map_t>
gpg_from_circuit(const Circuit& circ) {
  // circuit for conversion
  Circuit circ_flat(circ);
  unsigned n_qubits = circ_flat.n_qubits();
  unsigned n_bits = circ_flat.n_bits();
  // empty circuit
  Circuit empty_circ(n_qubits, n_bits);
  std::optional<std::string> name = circ_flat.get_name();
  if (name != std::nullopt) {
    empty_circ.set_name(name.value());
  }
  empty_circ.add_phase(circ_flat.get_phase());
  // measurement circuit
  Circuit measure_circ(n_qubits, n_bits);

  // flatten registers before process
  unit_map_t unit_map = circ_flat.flatten_registers();
  unit_map_t rev_unit_map;
  for (const auto& pair : unit_map) {
    rev_unit_map.insert({pair.second, pair.first});
  }

  std::vector<std::vector<PauliNode_ptr>> rotation_sets;
  std::vector<Command> commands = circ_flat.get_commands();
  Circuit cliff(n_qubits);
  // extract the final clifford and the measurement circuits
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::Measure: {
        measure_circ.add_op<UnitID>(OpType::Measure, cmd.get_args());
        break;
      }
      default: {
        if (optype == OpType::PauliExpBox ||
            optype == OpType::PauliExpPairBox ||
            optype == OpType::PauliExpCommutingSetBox)
          break;
        TKET_ASSERT(is_clifford_type(optype) && is_gate_type(optype));
        cliff.add_op<UnitID>(optype, cmd.get_args());
      }
    }
  }
  UnitaryRevTableau tab = circuit_to_unitary_rev_tableau(cliff);
  // convert the tableau into a set of nodes
  std::vector<PauliNode_ptr> rows = get_nodes_from_tableau(tab, n_qubits);
  // extract the Pauli exps
  for (const Command& cmd : commands) {
    OpType optype = cmd.get_op_ptr()->get_type();
    switch (optype) {
      case OpType::PauliExpBox: {
        const PauliExpBox& pbox =
            static_cast<const PauliExpBox&>(*cmd.get_op_ptr());
        const Expr phase = pbox.get_phase();
        const std::vector<Pauli> paulis = pbox.get_paulis();
        auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
        if (trivial) {
          empty_circ.add_phase(global_phase);
        } else {
          rotation_sets.push_back(
              {get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits)});
        }
        break;
      }
      case OpType::PauliExpPairBox: {
        const PauliExpPairBox& pbox =
            static_cast<const PauliExpPairBox&>(*cmd.get_op_ptr());
        const auto [paulis1, paulis2] = pbox.get_paulis_pair();
        const auto [phase1, phase2] = pbox.get_phase_pair();
        auto [trivial1, global_phase1] = is_trivial_pauliexp(paulis1, phase1);
        auto [trivial2, global_phase2] = is_trivial_pauliexp(paulis2, phase2);
        std::vector<PauliNode_ptr> rotation_set;
        if (trivial1) {
          empty_circ.add_phase(global_phase1);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis1, phase1, cmd.get_qubits(), n_qubits));
        }
        if (trivial2) {
          empty_circ.add_phase(global_phase2);
        } else {
          rotation_set.push_back(
              get_node_from_exp(paulis2, phase2, cmd.get_qubits(), n_qubits));
        }
        if (!rotation_set.empty()) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      case OpType::PauliExpCommutingSetBox: {
        const PauliExpCommutingSetBox& pbox =
            static_cast<const PauliExpCommutingSetBox&>(*cmd.get_op_ptr());
        const std::vector<SymPauliTensor> gadgets = pbox.get_pauli_gadgets();
        std::vector<PauliNode_ptr> rotation_set;
        for (const SymPauliTensor& pt : gadgets) {
          const std::vector<Pauli> paulis = pt.string;
          const Expr phase = pt.coeff;
          auto [trivial, global_phase] = is_trivial_pauliexp(paulis, phase);
          if (trivial) {
            empty_circ.add_phase(global_phase);
          } else {
            rotation_set.push_back(
                get_node_from_exp(paulis, phase, cmd.get_qubits(), n_qubits));
          }
        }
        if (rotation_set.size() > 0) {
          rotation_sets.push_back(rotation_set);
        }
        break;
      }
      default:
        break;
    }
  }
  return {empty_circ, rotation_sets, rows, measure_circ, rev_unit_map};
}

}  // namespace GreedyPauliSimp

}  // namespace Transforms

}  // namespace tket

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

#include "PauliOptimisation.hpp"

#include "Converters/Converters.hpp"
#include "Converters/PauliGadget.hpp"
#include "Decomposition.hpp"
#include "OptimisationPass.hpp"
#include "PauliGraph/PauliGraph.hpp"
#include "Transform.hpp"

namespace tket {

namespace Transforms {

Transform pairwise_pauli_gadgets(CXConfigType cx_config) {
  return Transform([=](Circuit &circ) {
    Expr t = circ.get_phase();
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType optype = op->get_type();
      if (is_boundary_c_type(optype)) continue;
      if (optype == OpType::Conditional)
        throw CircuitInvalidity(
            "Cannot currently do `pauli_gadgets` optimisation on a "
            "circuit with conditional gates");
      if (op->get_desc().is_box())
        throw CircuitInvalidity(
            "Cannot currently do `pauli_gadgets` optimisation on a "
            "circuit with boxes");
      if (optype == OpType::Measure || optype == OpType::Collapse) {
        VertexVec q_suc = circ.get_successors_of_type(v, EdgeType::Quantum);
        if (q_suc.size() != 1 ||
            !is_final_q_type(circ.get_OpType_from_Vertex(q_suc[0])))
          throw CircuitInvalidity(
              "Cannot currently do `pauli_gadgets` optimisation "
              "on a circuit with a measure in the middle of the "
              "circuit");
      }
    }
    Transform setup = decompose_multi_qubits_CX() >> decompose_ZX() >>
                      decompose_ZX_to_cliffords();
    setup.apply(circ);
    // We effectively commute non-Clifford rotations to the front of the circuit
    // This gives a sequence of just Pauli gadgets (gadget_circ), followed by
    // all of the Clifford operations (clifford_circ)
    std::vector<std::pair<QubitPauliTensor, Expr>> pauli_gadgets;
    // rx_pauli[i] specifies which Pauli gadget would be built by applying an Rx
    // rotation on qubit i and then pushing it through the Cliffords to the
    // front of the circuit. Likewise for rz_pauli with Rz rotations. Clifford
    // operations will update these and non-Clifford rotations will introduce
    // Pauli gadgets accordingly
    Circuit gadget_circ;
    Circuit clifford_circ;
    std::map<Qubit, QubitPauliTensor> rx_pauli;
    std::map<Qubit, QubitPauliTensor> rz_pauli;
    for (const Qubit &qb : circ.all_qubits()) {
      gadget_circ.add_qubit(qb);
      clifford_circ.add_qubit(qb);
      rx_pauli.insert({qb, QubitPauliTensor(qb, Pauli::X)});
      rz_pauli.insert({qb, QubitPauliTensor(qb, Pauli::Z)});
    }
    for (const Bit &cb : circ.all_bits()) {
      gadget_circ.add_bit(cb);
      clifford_circ.add_bit(cb);
    }
    // Identify Pauli Gadgets and build Clifford circuit
    for (const Command &c : circ) {
      const Op_ptr op_ptr = c.get_op_ptr();
      unit_vector_t args = c.get_args();
      OpType type = op_ptr->get_type();
      switch (type) {
        // Update rx_pauli and rz_pauli
        case OpType::S: {
          Qubit q(args[0]);
          rx_pauli[q] = i_ * rz_pauli[q] * rx_pauli[q];
          break;
        }
        case OpType::V: {
          Qubit q(args[0]);
          rz_pauli[q] = i_ * rx_pauli[q] * rz_pauli[q];
          break;
        }
        case OpType::Z: {
          Qubit q(args[0]);
          rx_pauli[q] = -1. * rx_pauli[q];
          break;
        }
        case OpType::X: {
          Qubit q(args[0]);
          rz_pauli[q] = -1. * rz_pauli[q];
          break;
        }
        case OpType::Sdg: {
          Qubit q(args[0]);
          rx_pauli[q] = -i_ * rz_pauli[q] * rx_pauli[q];
          break;
        }
        case OpType::Vdg: {
          Qubit q(args[0]);
          rz_pauli[q] = -i_ * rx_pauli[q] * rz_pauli[q];
          break;
        }
        case OpType::CX: {
          Qubit q_ctrl(args[0]);
          Qubit q_trgt(args[1]);
          rx_pauli[q_ctrl] = rx_pauli[q_ctrl] * rx_pauli[q_trgt];
          rz_pauli[q_trgt] = rz_pauli[q_ctrl] * rz_pauli[q_trgt];
          break;
        }
        // Introduce a Pauli gadget
        case OpType::Rz: {
          Qubit q(args[0]);
          Expr angle = (op_ptr)->get_params()[0];
          pauli_gadgets.push_back({rz_pauli[q], angle});
          break;
        }
        case OpType::Rx: {
          Qubit q(args[0]);
          Expr angle = (op_ptr)->get_params()[0];
          pauli_gadgets.push_back({rx_pauli[q], angle});
          break;
        }
        case OpType::Measure:
        case OpType::Collapse:
        case OpType::Reset:
          break;
        default: {
          std::string error_gate =
              "Cannot perform pairwise Pauli gadget optimisation "
              "using: " +
              op_ptr->get_name();
          throw NotValid(error_gate);
        }
      }
      // Add Clifford gates to the back of the circuit to recreate the final
      // combination at the outputs
      switch (type) {
        case OpType::Rz:
        case OpType::Rx: {
          break;
        }
        default: {
          clifford_circ.add_op(op_ptr, args);
        }
      }
    }
    // Synthesise pairs of Pauli Gadgets
    unsigned g = 0;
    while (g + 1 < pauli_gadgets.size()) {
      auto [pauli0, angle0] = pauli_gadgets[g];
      auto [pauli1, angle1] = pauli_gadgets[g + 1];
      append_pauli_gadget_pair(
          gadget_circ, pauli0, angle0, pauli1, angle1, cx_config);
      g += 2;
    }
    // As we synthesised Pauli gadgets 2 at a time, if there were an odd
    // number, we will have one left over, so add that one on its own
    if (g < pauli_gadgets.size()) {
      auto [pauli, angle] = pauli_gadgets[g];
      append_single_pauli_gadget(gadget_circ, pauli, angle, cx_config);
    }
    // Stitch gadget circuit and Clifford circuit together
    circ = gadget_circ >> clifford_circ;
    circ.add_phase(t);
    clifford_simp().apply(circ);
    return true;
  });
}

Transform synthesise_pauli_graph(
    PauliSynthStrat strat, CXConfigType cx_config) {
  return Transform([=](Circuit &circ) {
    Expr t = circ.get_phase();
    PauliGraph pg = circuit_to_pauli_graph(circ);
    switch (strat) {
      case PauliSynthStrat::Individual: {
        circ = pauli_graph_to_circuit_individually(pg, cx_config);
        break;
      }
      case PauliSynthStrat::Pairwise: {
        circ = pauli_graph_to_circuit_pairwise(pg, cx_config);
        break;
      }
      case PauliSynthStrat::Sets: {
        circ = pauli_graph_to_circuit_sets(pg, cx_config);
        break;
      }
      default:
        throw NotImplemented("Unknown Pauli Synthesis Strategy");
    }
    circ.add_phase(t);
    // always turn circuit into PauliGraph and back, so always return true
    return true;
  });
}

Transform special_UCC_synthesis(PauliSynthStrat strat, CXConfigType cx_config) {
  return Transform([=](Circuit &circ) {
    Transform synther = synthesise_pauli_graph(strat, cx_config);
    // make list so we don't run into unboxing vertex issues
    std::list<Vertex> circbox_verts;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::CircBox)
        circbox_verts.push_back(v);
    }
    for (Vertex v : circbox_verts) {
      const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      std::shared_ptr<const CircBox> box_ptr =
          std::dynamic_pointer_cast<const CircBox>(op);
      Circuit inner_circ = *(box_ptr->to_circuit());
      synther.apply(inner_circ);
      Subcircuit sub = {circ.get_in_edges(v), circ.get_all_out_edges(v), {v}};
      circ.substitute(inner_circ, sub);
    }
    return !circbox_verts
                .empty();  // always true if we have left Circuit formalism
  });
}

}  // namespace Transforms

}  // namespace tket

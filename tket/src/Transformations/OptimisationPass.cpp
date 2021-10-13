// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Circuit/CircUtils.hpp"
#include "Gate/GatePtr.hpp"
#include "Transform.hpp"

namespace tket {

Transform Transform::peephole_optimise_2q() {
  return (
      Transform::synthesise_tket() >> Transform::two_qubit_squash() >>
      Transform::hyper_clifford_squash() >> Transform::synthesise_tket());
}

Transform Transform::full_peephole_optimise(bool allow_swaps) {
  return (
      Transform::synthesise_tket() >> Transform::two_qubit_squash() >>
      Transform::clifford_simp(allow_swaps) >> Transform::synthesise_tket() >>
      Transform::three_qubit_squash() >>
      Transform::clifford_simp(allow_swaps) >> Transform::synthesise_tket());
}

Transform Transform::canonical_hyper_clifford_squash() {
  return Transform::optimise_via_PhaseGadget() >>
         Transform::two_qubit_squash() >> Transform::hyper_clifford_squash();
}

Transform Transform::hyper_clifford_squash() {
  return decompose_multi_qubits_CX() >> clifford_simp();
}

Transform Transform::clifford_simp(bool allow_swaps) {
  return decompose_cliffords_std() >> clifford_reduction(allow_swaps) >>
         decompose_multi_qubits_CX() >> singleq_clifford_sweep() >>
         squash_1qb_to_tk1();
}

Transform Transform::synthesise_tket() {
  Transform seq =
      Transform::commute_through_multis() >> Transform::remove_redundancies();
  Transform repeat = Transform::repeat(seq);
  Transform synth = Transform::decompose_multi_qubits_CX() >>
                    Transform::remove_redundancies() >> repeat >>
                    Transform::squash_1qb_to_tk1();
  Transform small_part = Transform::remove_redundancies() >> repeat >>
                         Transform::squash_1qb_to_tk1();
  Transform repeat_synth = Transform::repeat_with_metric(
      small_part, [](const Circuit &circ) { return circ.n_vertices(); });
  return synth >> repeat_synth;
}

static Transform CXs_from_phase_gadgets(CXConfigType cx_config) {
  return Transform([=](Circuit &circ) {
    bool success = false;
    VertexList bin;
    auto [i, end] = boost::vertices(circ.dag);
    for (auto next = i; i != end; i = next) {
      ++next;
      Vertex a = *i;
      Op_ptr op = circ.get_Op_ptr_from_Vertex(a);
      if (op->get_type() == OpType::PhaseGadget) {
        unsigned n_qubits = op->n_qubits();
        Circuit replacement =
            phase_gadget(n_qubits, op->get_params()[0], cx_config);
        Subcircuit sub = {
            {circ.get_in_edges(a)}, {circ.get_all_out_edges(a)}, {a}};
        circ.substitute(replacement, sub, Circuit::VertexDeletion::Yes);
        success = true;
      }
    }
    return success;
  });
}

Transform Transform::optimise_via_PhaseGadget(CXConfigType cx_config) {
  return rebase_tket() >> decompose_PhaseGadgets() >> smash_CX_PhaseGadgets() >>
         align_PhaseGadgets() >> CXs_from_phase_gadgets(cx_config) >>
         synthesise_tket();
}

Transform Transform::synthesise_OQC() {
  return Transform([](Circuit &circ) {
    Transform rep_zx = squash_1qb_to_pqp(OpType::Rx, OpType::Rz) >>
                       commute_through_multis() >> remove_redundancies();
    Transform seq = decompose_multi_qubits_CX() >> decompose_CX_to_ECR() >>
                    decompose_ZX() >> repeat(rep_zx) >> rebase_OQC() >>
                    commute_through_multis() >> remove_redundancies();
    return seq.apply(circ);
  });
}

/* Returns a Circuit with only HQS allowed Ops (Rz, PhasedX, ZZMax) */
Transform Transform::synthesise_HQS() {
  return Transform([](Circuit &circ) {
    Transform single_loop =
        remove_redundancies() >> commute_through_multis() >> reduce_XZ_chains();
    Transform hqs_loop = remove_redundancies() >> commute_and_combine_HQS2() >>
                         reduce_XZ_chains();
    Transform main_seq = decompose_multi_qubits_CX() >> clifford_simp() >>
                         decompose_ZX() >> repeat(single_loop) >>
                         decompose_CX_to_HQS2() >> repeat(hqs_loop) >>
                         decompose_ZX_to_HQS1();
    return main_seq.apply(circ);
  });
}

// TODO: Make the XXPhase gates combine
Transform Transform::synthesise_UMD() {
  return Transform([](Circuit &circ) {
    bool success = (Transform::synthesise_tket() >> Transform::decompose_ZX() >>
                    Transform::decompose_MolmerSorensen() >>
                    Transform::squash_1qb_to_tk1())
                       .apply(circ);
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      const Op_ptr op_ptr = circ.get_Op_ptr_from_Vertex(v);
      OpType type = op_ptr->get_type();
      if (type == OpType::tk1) {
        std::vector<Expr> tk1_angles = as_gate_ptr(op_ptr)->get_tk1_angles();
        Circuit in_circ = Transform::tk1_to_PhasedXRz(
            tk1_angles[0], tk1_angles[1], tk1_angles[2]);
        Subcircuit sub = {
            {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
        bin.push_back(v);
        circ.substitute(in_circ, sub, Circuit::VertexDeletion::No);
        circ.add_phase(tk1_angles[3]);
        success = true;
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

// Transform Transform::route_for_arc(const Architecture &arc, const bool
// decompose, const bool directed, const bool placement) {
//     return Transform([&](Circuit &circ) {
//         Routing route(circ,arc,decompose,directed,placement);
//         route.solve();
//         circ = route.returnCircuit();
//     });
// }

}  // namespace tket

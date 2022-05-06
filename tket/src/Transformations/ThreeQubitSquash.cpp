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

#include "ThreeQubitSquash.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <set>
#include <string>

#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Circuit/QInteraction.hpp"
#include "Circuit/ThreeQubitConversion.hpp"
#include "OpType/EdgeType.hpp"
#include "OpType/OpType.hpp"
#include "OptimisationPass.hpp"
#include "Transform.hpp"
#include "Utils/Assert.hpp"
#include "Utils/GraphHeaders.hpp"

namespace tket {

namespace Transforms {
// Candidate substitution for a 2- or 3-qubit circuit.
static Circuit candidate_sub(const Circuit &circ) {
  unsigned n_qb = circ.n_qubits();
  if (n_qb == 2) {
    Circuit repl = two_qubit_canonical(get_matrix_from_2qb_circ(circ));
    clifford_simp(false).apply(repl);
    return repl;
  } else {
    TKET_ASSERT(n_qb == 3);
    Circuit repl = three_qubit_synthesis(get_3q_unitary(circ));
    clifford_simp(false).apply(repl);
    return repl;
  }
}

Transform three_qubit_squash() {
  return Transform([](Circuit &circ) {
    bool changed = false;

    // Step through the vertices in topological order.
    QISystem Is(circ, candidate_sub);  // set of "live" interactions
    for (const Vertex &v : circ.vertices_in_order()) {
      const EdgeVec v_q_ins = circ.get_in_edges_of_type(v, EdgeType::Quantum);
      const EdgeVec v_q_outs = circ.get_out_edges_of_type(v, EdgeType::Quantum);
      unsigned n_q_ins = v_q_ins.size();
      unsigned n_q_outs = v_q_outs.size();

      // If v has no quantum wires, ignore it and move on.
      if (n_q_ins == 0 && n_q_outs == 0) continue;

      // If v is initial, create an interaction from its out-edge, and move on.
      if (n_q_ins == 0) {
        TKET_ASSERT(n_q_outs == 1);
        Is.create_new_interaction_from_edge(v_q_outs[0]);
        continue;
      }

      // If v is final, ignore it and move on.
      if (n_q_outs == 0) continue;

      // It's an internal operation with >0 quantum wires.
      TKET_ASSERT(n_q_ins == n_q_outs);

      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType optype = op->get_type();

      // If there are any incoming classical wires, or if this is a Barrier or
      // Reset or Collapse operation, or if the operation contains symbols,
      // close all existing interactions meeting v and create new ones, then
      // move on.
      if (!circ.get_in_edges_of_type(v, EdgeType::Classical).empty() ||
          !circ.get_in_edges_of_type(v, EdgeType::Boolean).empty() ||
          optype == OpType::Barrier || optype == OpType::Reset ||
          optype == OpType::Collapse || !op->free_symbols().empty()) {
        changed |= Is.close_interactions_feeding_vertex(v);
        continue;
      }

      // The circuit should contain only 1-qubit and CX gates.
      if ((n_q_ins == 2 && optype != OpType::CX) || (n_q_ins > 2)) {
        throw std::invalid_argument(
            "Three-qubit squash requires circuits with 1q and CX gates only");
      }

      // Absorb v into existing interactions, closing or merging as necessary.
      bool done_with_v = false;
      while (!done_with_v) {
        std::vector<int> v_Is = Is.interactions_feeding_vertex(v);
        unsigned total_n_wires = Is.total_n_wires(v_Is);
        if (total_n_wires <= 3) {
          Is.combine_and_append(v_Is, v);
          done_with_v = true;
        } else {
          // Close one of the interactions meeting v.
          int i = Is.largest_interaction(v_Is);
          changed |= Is.close_interaction_and_spawn(i);
        }
      }
    }

    // Close all remaining interactions.
    changed |= Is.close_all_interactions();

    // Delete removed vertices.
    Is.destroy_bin();

    return changed;
  });
}

}  // namespace Transforms

}  // namespace tket

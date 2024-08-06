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

#include "tket/Transformations/CliffordResynthesis.hpp"

#include <algorithm>
#include <vector>

#include "tket/Circuit/Circuit.hpp"
#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/OptimisationPass.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

static Circuit resynthesised_circuit(
    const Circuit &circ,
    std::optional<std::function<Circuit(const Circuit &)>> transform = nullptr,
    bool allow_swaps = true) {
  if (transform.has_value()) {
    return (*transform)(circ);
  } else {
    Circuit circ1(circ);
    // Ensure all recognized Clifford gates are Clifford types:
    decompose_multi_qubits_CX().apply(circ1);
    decompose_cliffords_std().apply(circ1);
    remove_redundancies().apply(circ1);

    // Convert to tableau and back:
    const UnitaryTableau tab = circuit_to_unitary_tableau(circ1);
    Circuit newcirc = unitary_tableau_to_circuit(tab);

    // Simplify:
    clifford_simp(allow_swaps).apply(newcirc);

    return newcirc;
  }
}

// Returns the two-qubit-gate reduction resulting from resynthesising the given
// (convex, connected) vertex set using the given transform.
static int n_2q_reduction(
    const Circuit &circ, const VertexSet &verts,
    std::optional<std::function<Circuit(const Circuit &)>> transform,
    bool allow_swaps) {
  Subcircuit subcircuit = circ.make_subcircuit(verts);
  const Circuit subc = circ.subcircuit(subcircuit);
  Circuit newsubc = resynthesised_circuit(subc, transform, allow_swaps);
  return int(subc.count_n_qubit_gates(2)) - int(newsubc.count_n_qubit_gates(2));
}

static bool resynthesise_cliffords(
    Circuit &circ,
    std::optional<std::function<Circuit(const Circuit &)>> transform = nullptr,
    bool allow_swaps = true) {
  bool changed = false;
  while (true) {
    std::vector<VertexSet> subcircuits =
        circ.get_subcircuits([](Op_ptr op) { return op->is_clifford(); });
    if (subcircuits.empty()) {
      break;
    }
    // Pick the one that reduces the 2q gate count the most.
    VertexSet verts = *std::max_element(
        subcircuits.begin(), subcircuits.end(),
        [&circ, &transform, &allow_swaps](
            const VertexSet &A, const VertexSet &B) {
          return n_2q_reduction(circ, A, transform, allow_swaps) <
                 n_2q_reduction(circ, B, transform, allow_swaps);
        });
    Subcircuit subcircuit = circ.make_subcircuit(verts);
    const Circuit subc = circ.subcircuit(subcircuit);
    Circuit newsubc = resynthesised_circuit(subc, transform, allow_swaps);
    if (newsubc.count_n_qubit_gates(2) < subc.count_n_qubit_gates(2)) {
      circ.substitute(newsubc, subcircuit);
      changed = true;
    } else {
      break;
    }
  }
  return changed;
}

Transform clifford_resynthesis(
    std::optional<std::function<Circuit(const Circuit &)>> transform,
    bool allow_swaps) {
  return Transform([transform, allow_swaps](Circuit &circ) {
    return resynthesise_cliffords(circ, transform, allow_swaps);
  });
}

}  // namespace Transforms

}  // namespace tket

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

#include "Transformations/Decomposition.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Clifford/UnitaryTableau.hpp"
#include "tket/Converters/Converters.hpp"
#include "tket/Transformations/OptimisationPass.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

static bool resynthesise_cliffords(Circuit &circ, bool allow_swaps = true) {
  bool changed = false;
  for (const VertexSet &verts :
       circ.get_subcircuits([](Op_ptr op) { return op->is_clifford(); })) {
    Subcircuit subcircuit = circ.make_subcircuit(verts);
    Circuit subc = circ.subcircuit(subcircuit);
    decompose_cliffords_std().apply(
        subc);  // works around https://github.com/CQCL/tket/issues/1268
    const UnitaryTableau tab = circuit_to_unitary_tableau(subc);
    Circuit newsubc = unitary_tableau_to_circuit(tab);
    clifford_simp(allow_swaps).apply(newsubc);
    circ.substitute(newsubc, subcircuit);
    changed = true;
  }
  return changed;
}

Transform clifford_resynthesis(bool allow_swaps) {
  return Transform([allow_swaps](Circuit &circ) {
    return resynthesise_cliffords(circ, allow_swaps);
  });
}

}  // namespace Transforms

}  // namespace tket

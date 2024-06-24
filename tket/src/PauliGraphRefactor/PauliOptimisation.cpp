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

#include "tket/PauliGraphRefactor/PauliOptimisation.hpp"

#include "tket/PauliGraphRefactor/Converters.hpp"

namespace tket {

using namespace pg;

namespace Transforms {

Transform synthesise_pauli_graph3(
    PauliSynthStrat3 strat, CXConfigType cx_config) {
  return Transform([=](Circuit &circ) {
    std::optional<std::string> name = circ.get_name();
    PauliGraph pg = circuit_to_pauli_graph3(circ);
    switch (strat) {
      case PauliSynthStrat3::Individual: {
        circ =
            pauli_graph3_to_pauli_exp_box_circuit_individually(pg, cx_config);
        break;
      }
      case PauliSynthStrat3::Pairwise: {
        circ = pauli_graph3_to_pauli_exp_box_circuit_pairwise(pg, cx_config);
        break;
      }
      case PauliSynthStrat3::Sets: {
        circ = pauli_graph3_to_pauli_exp_box_circuit_sets(pg, cx_config);
        break;
      }
      case PauliSynthStrat3::ExtIndividual: {
        circ = pauli_graph3_to_circuit_individual(pg, cx_config);
        break;
      }
      case PauliSynthStrat3::ExtSets: {
        circ = pauli_graph3_to_circuit_sets(pg, cx_config);
        break;
      }
      default:
        TKET_ASSERT(!"Unknown Pauli Synthesis Strategy");
    }
    if (name) circ.set_name(*name);
    // always turn circuit into PauliGraph and back, so always return true
    return true;
  });
}

}  // namespace Transforms

}  // namespace tket

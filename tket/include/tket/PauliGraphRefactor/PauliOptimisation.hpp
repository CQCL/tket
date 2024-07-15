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

#pragma once

#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

/* Dictates whether synthesis of a PauliGraph should
   be done on the Paulis individually, making use of the pairwise
   interactions or collecting into mutually commuting sets. Divided into legacy
   options (only works for unitary circuits with end-of-circuit measurements) or
   extended PauliGraph options. */
enum class PauliSynthStrat3 {
  Individual,
  Pairwise,
  Sets,
  ExtIndividual,
  ExtSets
};

NLOHMANN_JSON_SERIALIZE_ENUM(
    PauliSynthStrat3, {{PauliSynthStrat3::Individual, "Individual"},
                       {PauliSynthStrat3::Pairwise, "Pairwise"},
                       {PauliSynthStrat3::Sets, "Sets"},
                       {PauliSynthStrat3::ExtIndividual, "ExtIndividual"},
                       {PauliSynthStrat3::ExtSets, "ExtSets"}});

// always returns true, as it leaves Circuit data structure
Transform synthesise_pauli_graph3(
    PauliSynthStrat3 strat = PauliSynthStrat3::Sets,
    CXConfigType cx_config = CXConfigType::Snake);

}  // namespace Transforms

}  // namespace tket

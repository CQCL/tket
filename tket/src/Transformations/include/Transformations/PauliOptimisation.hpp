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

#pragma once

#include "Transform.hpp"

namespace tket {

namespace Transforms {

/* Dictates whether synthesis of a PauliGraph should
    be done on the Paulis individually, making use of the pairwise
    interactions or collecting into mutually commuting sets. */
enum class PauliSynthStrat { Individual, Pairwise, Sets };

NLOHMANN_JSON_SERIALIZE_ENUM(
    PauliSynthStrat, {{PauliSynthStrat::Individual, "Individual"},
                      {PauliSynthStrat::Pairwise, "Pairwise"},
                      {PauliSynthStrat::Sets, "Sets"}});

Transform pairwise_pauli_gadgets(CXConfigType cx_config = CXConfigType::Snake);

// always returns true, as it leaves Circuit data structure
Transform synthesise_pauli_graph(
    PauliSynthStrat strat = PauliSynthStrat::Sets,
    CXConfigType cx_config = CXConfigType::Snake);

// Assumes incoming circuit is composed of `CircBox`es with
// `PauliExpBox`es inside
Transform special_UCC_synthesis(
    PauliSynthStrat strat = PauliSynthStrat::Sets,
    CXConfigType cx_config = CXConfigType::Snake);

}  // namespace Transforms

}  // namespace tket

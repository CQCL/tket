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

#include <catch2/catch.hpp>
#include <numeric>
#include <optional>

#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "CircuitsForTesting.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpType.hpp"
#include "Ops/MetaOp.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/CliffordOptimisation.hpp"
#include "Transformations/Combinator.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Replacement.hpp"
#include "Transformations/Transform.hpp"
#include "testutil.hpp"

/* This test file covers decomposition, basic optimisation and synthesis passes.
It does not cover Rebasing, Clifford optimisation, Phase Gadgets,
Multi-controlled Decomp, CZ Optimisation, PauliString optimisation, getting
matrices from Circuits etc*/

// TODO: Split this up more

namespace tket {
namespace test_Synthesis {

SCENARIO("Testing globalise_PhasedX") {
  GIVEN("A more complex NPhasedX circuit on 4qb circuit") {
    Circuit circ(3);
    // circ.add_op<unsigned>(OpType::NPhasedX, {0.2, 0.54}, {0, 1});
    circ.add_op<unsigned>(OpType::NPhasedX, {0.53, 0.23}, {0, 1});
    // circ.add_op<unsigned>(OpType::PhasedX, {0.3, 0.2}, {1});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    // circ.add_op<unsigned>(OpType::CX, {0, 2});
    // circ.add_op<unsigned>(OpType::H, {1});
    // circ.add_op<unsigned>(OpType::H, {3});
    circ.add_op<unsigned>(OpType::NPhasedX, {0.53, 0.23}, {0, 1, 2});
    Circuit tmp_circ(circ);
    Transforms::globalise_PhasedX().apply(tmp_circ);
  }
}

}  // namespace test_Synthesis
}  // namespace tket

// Copyright 2019-2023 Cambridge Quantum Computing
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

#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <random>

#include "../testutil.hpp"
#include "Circuit/Boxes.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Circuit.hpp"
#include "Circuit/StatePermutation.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Gate/Rotation.hpp"
#include "Simulation/CircuitSimulator.hpp"

namespace tket {
namespace test_DiagonalBox {


SCENARIO("Test StatePermutationBox") {
    GIVEN("2-q permutation"){
        std::map<std::vector<bool>, std::vector<bool>> permutation;

        permutation[{0, 0}] = {0, 0};
        permutation[{0, 1}] = {1, 1};
        permutation[{1, 0}] = {0, 1};
        permutation[{1, 1}] = {1, 0};

        StatePermutationBox box(permutation);
        Circuit circ = *box.to_circuit();
        const auto matrix = tket_sim::get_unitary(circ);

        ToffoliBox tb(2, permutation);
        Circuit circ2 = *tb.to_circuit();
        const auto matrix2 = tket_sim::get_unitary(circ2);

        REQUIRE((matrix - matrix2).cwiseAbs().sum() < ERR_EPS);
    }
    GIVEN("3-q permutation"){
        std::map<std::vector<bool>, std::vector<bool>> permutation;
        permutation[{0, 0, 0}] = {1, 0, 0};
        permutation[{0, 0, 1}] = {0, 0, 1};
        permutation[{0, 1, 0}] = {1, 0, 1};
        permutation[{0, 1, 1}] = {0, 1, 0};
        permutation[{1, 0, 0}] = {0, 0, 0};
        permutation[{1, 0, 1}] = {0, 1, 1};
        permutation[{1, 1, 0}] = {1, 1, 1};
        permutation[{1, 1, 1}] = {1, 1, 0};

        StatePermutationBox box(permutation);
        Circuit circ = *box.to_circuit();
        const auto matrix = tket_sim::get_unitary(circ);

        ToffoliBox tb(3, permutation);
        Circuit circ2 = *tb.to_circuit();
        const auto matrix2 = tket_sim::get_unitary(circ2);
        REQUIRE((matrix - matrix2).cwiseAbs().sum() < ERR_EPS);
    }
}

}  // namespace test_DiagonalBox
}  // namespace tket

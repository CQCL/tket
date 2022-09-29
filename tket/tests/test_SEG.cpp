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

#include <catch2/catch_test_macros.hpp>
#include <numeric>
#include <optional>

#include "Characterisation/DeviceCharacterisation.hpp"
#include "Circuit/Circuit.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "Mapping/Verification.hpp"
#include "OpType/OpType.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassGenerators.hpp"
#include "Predicates/Predicates.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/HelperFunctions.hpp"
#include "testutil.hpp"

namespace tket {
SCENARIO("USER SEG") {
  GIVEN("MINIMMUM VIABLE") {
    Circuit circ(25);
    add_2qb_gates(
        circ, OpType::CX,
        {{2, 1},
         {3, 7},
         {0, 3},
         {6, 9},
         {7, 15},
         {16, 6},
         {18, 12},
         {7, 19},
         {4, 21},
         {18, 4},
         {23, 11},
         {17, 24},
         {8, 13}});
    circ.add_barrier({0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12,
                      13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
    add_2qb_gates(circ, OpType::CX, {{2, 1}, {23, 19}, {23, 11}});

    std::vector<std::pair<unsigned, unsigned>> edges = {
        {0, 1},   {0, 5},   {0, 6},   {1, 0},   {1, 2},   {1, 5},   {1, 6},
        {1, 7},   {2, 1},   {2, 3},   {2, 6},   {2, 7},   {2, 8},   {3, 2},
        {3, 4},   {3, 7},   {3, 8},   {3, 9},   {4, 3},   {4, 8},   {4, 9},
        {5, 0},   {5, 1},   {5, 6},   {5, 10},  {5, 11},  {6, 0},   {6, 1},
        {6, 2},   {6, 5},   {6, 7},   {6, 10},  {6, 11},  {6, 12},  {7, 1},
        {7, 2},   {7, 3},   {7, 6},   {7, 8},   {7, 11},  {7, 12},  {7, 13},
        {8, 2},   {8, 3},   {8, 4},   {8, 7},   {8, 9},   {8, 12},  {8, 13},
        {8, 14},  {9, 3},   {9, 4},   {9, 8},   {9, 13},  {9, 14},  {10, 5},
        {10, 6},  {10, 11}, {10, 15}, {10, 16}, {11, 5},  {11, 6},  {11, 7},
        {11, 10}, {11, 12}, {11, 15}, {11, 16}, {11, 17}, {12, 6},  {12, 7},
        {12, 8},  {12, 11}, {12, 13}, {12, 16}, {12, 17}, {12, 18}, {13, 7},
        {13, 8},  {13, 9},  {13, 12}, {13, 14}, {13, 17}, {13, 18}, {13, 19},
        {14, 8},  {14, 9},  {14, 13}, {14, 18}, {14, 19}, {15, 10}, {15, 11},
        {15, 16}, {15, 20}, {15, 21}, {16, 10}, {16, 11}, {16, 12}, {16, 15},
        {16, 17}, {16, 20}, {16, 21}, {16, 22}, {17, 11}, {17, 12}, {17, 13},
        {17, 16}, {17, 18}, {17, 21}, {17, 22}, {17, 23}, {18, 12}, {18, 13},
        {18, 14}, {18, 17}, {18, 19}, {18, 22}, {18, 23}, {18, 24}, {19, 13},
        {19, 14}, {19, 18}, {19, 23}, {19, 24}, {20, 15}, {20, 16}, {20, 21},
        {21, 15}, {21, 16}, {21, 17}, {21, 20}, {21, 22}, {22, 16}, {22, 17},
        {22, 18}, {22, 21}, {22, 23}, {23, 17}, {23, 18}, {23, 19}, {23, 22},
        {23, 24}, {24, 18}, {24, 19}, {24, 23},
    };
    Architecture arc(edges);
    PassPtr r_p = gen_routing_pass(
        arc, {std::make_shared<LexiLabellingMethod>(),
              std::make_shared<LexiRouteRoutingMethod>()});
    CompilationUnit cu(circ);
    r_p->apply(cu);
    std::cout << cu.get_circ_ref() << std::endl;
  }
}

}  // namespace tket
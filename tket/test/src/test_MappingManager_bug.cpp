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

#include <catch2/catch_test_macros.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "ArchAwareSynth/SteinerForest.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/Circuit.hpp"
#include "Converters/PhasePoly.hpp"
#include "Mapping/LexiLabelling.hpp"
#include "Mapping/LexiRoute.hpp"
#include "Mapping/MappingManager.hpp"
#include "OpType/OpType.hpp"
#include "Placement/Placement.hpp"
#include "Predicates/CompilationUnit.hpp"
#include "Predicates/CompilerPass.hpp"
#include "Predicates/PassLibrary.hpp"
#include "Predicates/Predicates.hpp"
#include "Transformations/BasicOptimisation.hpp"
#include "Transformations/ContextualReduction.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/PauliOptimisation.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/ThreeQubitSquash.hpp"
#include "Transformations/Transform.hpp"
#include "Utils/Json.hpp"


namespace tket {

SCENARIO("Segfault Github #777") {
    std::ifstream arch_file("ibm_montreal.json");
    nlohmann::json j_arch = nlohmann::json::parse(arch_file);
    auto arch = j_arch.get<Architecture>();
    std::ifstream circ_file("bug777_circuit.json");
    nlohmann::json j_circ = nlohmann::json::parse(circ_file);
    auto circ = j_circ.get<Circuit>();
    std::map<Qubit, Node> p_map = {
        {Node(0), Node(5)},
        {Node(1), Node(8)},
        {Node(2), Node("unplaced", 0)},
        {Node(3), Node(16)},
        {Node(4), Node(3)},
        {Node(5), Node("unplaced", 1)},
        {Node(6), Node("unplaced", 2)},
        {Node(7), Node("unplaced", 3)},
        {Node(8), Node("unplaced", 4)},
        {Node(9), Node("unplaced", 5)},
        {Node(10), Node("unplaced", 6)},
        {Node(11), Node(25)},
        {Node(12), Node("unplaced", 7)},
        {Node(13), Node(14)},
        {Node(14), Node("unplaced", 8)},
        {Node(15), Node(19)},
        {Node(16), Node(24)},
        {Node(17), Node("unplaced", 9)},
        {Node(18), Node("unplaced", 10)},
        {Node(19), Node(2)},
        {Node(20), Node(1)},
        {Node(21), Node(22)},
        {Node(22), Node(11)},
        {Node(23), Node("unplaced", 11)},
        {Node(24), Node("unplaced", 12)},
        {Node(25), Node("unplaced", 13)},
        {Node(26), Node("unplaced", 14)}};
    MappingManager mm(std::make_shared<Architecture>(arch));
    std::shared_ptr<unit_bimaps_t> maps = std::make_shared<unit_bimaps_t>();
    for (auto it : p_map){
      maps->initial.insert({it.first, it.second});
      maps->final.insert({it.first, it.second});
    }
    std::vector<RoutingMethodPtr> config = {
        std::make_shared<LexiLabellingMethod>(),
        std::make_shared<LexiRouteRoutingMethod>()};
    REQUIRE(mm.route_circuit_with_maps(circ, config, maps));
}
}  // namespace tket

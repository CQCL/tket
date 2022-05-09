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
#include <random>

#include "WeightSubgrMono/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("Test search node string functions") {
  PossibleAssignments possible_assignments;
  possible_assignments[0] = {0, 1};
  possible_assignments[3] = {2};
  NodesRawData nodes_raw_data(possible_assignments);
  auto& node_data = nodes_raw_data.nodes_data.at(0);

  node_data.new_assignments.emplace_back(0, 0);
  CHECK(
      node_data.str() ==
      "Has 2 ass.: [ 3:2 0:0 ];  sc.prod 0; p-edge weight 0");

  CHECK(nodes_raw_data.domains_data.at(3).str() == "\n  lev=0, Dom: [ 2 ]\n");

  node_data.nogood = true;
  CHECK(
      node_data.str() ==
      "##NOGOOD!## Has 2 ass.: [ 3:2 0:0 ];  sc.prod 0; p-edge weight 0");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

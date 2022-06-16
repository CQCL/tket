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
#include <stdexcept>

#include "WeightSubgrMono/Searching/NodesRawData.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("Test search node string functions") {
  DomainInitialiser::InitialDomains initial_domains(4);
  initial_domains[0] = {0, 1};
  initial_domains[3] = {2};
  CHECK_THROWS_AS(NodesRawData(initial_domains), std::runtime_error);
  initial_domains[1] = {17};
  initial_domains[2] = {77, 88};

  NodesRawData nodes_raw_data(initial_domains);
  auto& node_data = nodes_raw_data.nodes_data[0];

  node_data.new_assignments.emplace_back(0, 0);
  CHECK(
      node_data.str() ==
      "Has 3 ass.: [ 1:17 3:2 0:0 ];  sc.prod 0; p-edge weight 0");

  CHECK(
      nodes_raw_data.domains_data.at(3).str() ==
      "\n  node_index=0, Dom: [ 2 ]\n");

  node_data.nogood = true;
  CHECK(
      node_data.str() ==
      "##NOGOOD!## Has 3 ass.: [ 1:17 3:2 0:0 ];  sc.prod 0; p-edge weight 0");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

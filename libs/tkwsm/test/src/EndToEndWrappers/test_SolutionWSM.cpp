// Copyright Quantinuum
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
#include <tkwsm/EndToEndWrappers/SolutionWSM.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

SCENARIO("SolutionWSM : errors with empty assignments but nonzero weights") {
  SolutionWSM solution;
  const GraphEdgeWeights pattern_edges_and_weights;
  const GraphEdgeWeights target_edges_and_weights;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) == "");

  solution.scalar_product = 10;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "empty assignments, but sc.prod=10, total p.edge.weights=0");

  solution.scalar_product = 0;
  solution.total_p_edges_weight = 20;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "empty assignments, but sc.prod=0, total p.edge.weights=20");
}

SCENARIO(
    "SolutionWSM : errors with nonempty but invalid/mismatching assignments") {
  GraphEdgeWeights pattern_edges_and_weights;
  pattern_edges_and_weights[get_edge(0, 1)] = 3;

  GraphEdgeWeights target_edges_and_weights;
  target_edges_and_weights[get_edge(0, 1)] = 5;
  target_edges_and_weights[get_edge(1, 2)] = 7;

  SolutionWSM solution;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) == "");

  solution.assignments.emplace_back(0, 1);
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "\nP-edge (0,1) has unassigned vertices\n"
      "Weights mismatch: scalar products 0,0; total p-edge weights 3,0\n"
      "Number of used p vertices mismatch: 2,1");

  solution.assignments.emplace_back(1, 2);
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "\nWeights mismatch: scalar products 21,0; total p-edge weights 3,0");

  solution.scalar_product = 21;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "\nWeights mismatch: scalar products 21,21; total p-edge weights 3,0");

  solution.total_p_edges_weight = 3;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) == "");

  solution.assignments.back().second = 999;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "\nP-edge [0,1] maps to nonexistent target edge [1,999]"
      "\nWeights mismatch: scalar products 0,21; total p-edge weights 3,3");

  solution.assignments.back().second = 1;
  CHECK(
      solution.get_errors(
          pattern_edges_and_weights, target_edges_and_weights) ==
      "\nDuplicate value 1 seen, when adding 1->1"
      "\nSizes mismatch: 2,2,1"
      "\nP vertices 0,1 both map to 1"
      "\nWeights mismatch: scalar products 0,21; total p-edge weights 3,3");
}

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

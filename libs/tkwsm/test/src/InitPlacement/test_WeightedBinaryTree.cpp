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

#include <catch2/catch_test_macros.hpp>
#include <tkwsm/Common/GeneralUtils.hpp>

#include "TestUtilsIQP.hpp"
#include "WeightedBinaryTree.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {
namespace InitialPlacement {
namespace tests {

// clang-format off
/* The fixed weighted binary tree will be (W means weight):

             1
           /   \
          /     \
         W1     W2
        /         \
       /           \
      2             3
     / \           / \
    W5  W7       W17 W30
   /     \       /     \
  4       5     6       7

*/
// clang-format on
SCENARIO("Paths and token swaps on small fixed binary tree") {
  const std::vector<WeightWSM> weights{0, 0, 1, 2, 5, 7, 17, 30};
  WeightedBinaryTree tree(weights, 3);

  // Manually checked
  CHECK(
      str(tree.get_graph_data()) ==
      "6 edges with weights: [  (1,2: 1),"
      "  (1,3: 2),  (2,4: 5),  (2,5: 7),  (3,6: 17),  (3,7: 30), ]\n"
      "7 vertices: {1 2 3 4 5 6 7 }\n");

  const std::vector<std::pair<VertexWSM, VertexWSM>> placement{
      {0, 1}, {1, 2}, {2, 4}, {3, 5}, {4, 6}, {5, 7}};

  tree.initialise_with_qubit_placement(placement);
  const auto& placement_map = tree.get_current_placement();
  REQUIRE(placement_map.size() == placement.size());

  const auto& tokens_map = tree.get_current_tokens();
  for (const auto& entry : placement) {
    REQUIRE(placement_map.at(entry.first) == entry.second);
    REQUIRE(tokens_map.at(entry.second) == entry.first);
  }

  const std::vector<std::pair<VertexWSM, VertexWSM>> gates{
      {2, 3}, {1, 5}, {2, 4}, {3, 0}, {4, 2}, {1, 2}};

  // clang-format off
  // Manually checked
  CHECK(do_token_swaps_and_check_placements(gates, tree) ==
    "\n"
    "TOKEN Swap (2,3) between vertices 4 5; cost 22\n"
    "Path vertices: [ 4 2 5 ]\n"
    "Edge weights: [ 5 7 ] (total weight 12)\n"
    "NOW, placement: { 0->1 1->4 2->2 3->5 4->6 5->7 }\n"
    "\n"
    "TOKEN Swap (1,5) between vertices 4 7; cost 54\n"
    "Path vertices: [ 4 2 1 3 7 ]\n"
    "Edge weights: [ 5 1 2 30 ] (total weight 38)\n"
    "NOW, placement: { 0->2 1->3 2->4 3->5 4->6 5->7 }\n"
    "\n"
    "TOKEN Swap (2,4) between vertices 4 6; cost 41\n"
    "Path vertices: [ 4 2 1 3 6 ]\n"
    "Edge weights: [ 5 1 2 17 ] (total weight 25)\n"
    "NOW, placement: { 0->4 1->1 2->3 3->5 4->6 5->7 }\n"
    "\n"
    "TOKEN Swap (3,0) between vertices 5 4; cost 22\n"
    "Path vertices: [ 5 2 4 ]\n"
    "Edge weights: [ 7 5 ] (total weight 12)\n"
    "NOW, placement: { 0->2 1->1 2->3 3->5 4->6 5->7 }\n"
    "\n"
    "TOKEN Swap (4,2) between vertices 6 3; cost 17\n"
    "Path vertices: [ 6 3 ]\n"
    "Edge weights: [ 17 ] (total weight 17)\n"
    "NOW, placement: { 0->2 1->1 2->3 3->5 4->6 5->7 }\n"
    "\n"
    "TOKEN Swap (1,2) between vertices 1 3; cost 2\n"
    "Path vertices: [ 1 3 ]\n"
    "Edge weights: [ 2 ] (total weight 2)\n"
    "NOW, placement: { 0->2 1->1 2->3 3->5 4->6 5->7 }\n"
  );
  // clang-format on
}

}  // namespace tests
}  // namespace InitialPlacement
}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

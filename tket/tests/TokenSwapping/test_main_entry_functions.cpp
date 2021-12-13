#include <catch2/catch.hpp>
#include <sstream>

#include "TokenSwapping/RNG.hpp"
#include "TokenSwapping/main_entry_functions.hpp"

;
using std::vector;

// Detailed algorithmic checks with quantitative benchmarks
// are done elsewhere, so this is really just checking conversion.

namespace tket {
namespace tsa_internal {
namespace tests {

SCENARIO("main entry function for TSA") {
  // Will summarise relevant data, so that we can see any changes.
  std::stringstream problem_ss;

  const SquareGrid arch(3, 4, 2);
  const auto nodes = arch.get_all_nodes_vec();
  const auto edges = arch.get_all_edges_vec();
  problem_ss << nodes.size() << " nodes; " << edges.size() << " edges.";

  // The value is the set of all neighbouring nodes.
  std::map<Node, std::set<Node>> allowed_edges_map;
  for (auto [n1, n2] : edges) {
    REQUIRE(n1 != n2);
    allowed_edges_map[n1].insert(n2);
    allowed_edges_map[n2].insert(n1);
  }

  // Key: a node Value: its original position in "nodes"
  std::map<Node, size_t> original_vertex_indices;
  for (size_t ii = 0; ii < nodes.size(); ++ii) {
    original_vertex_indices[nodes[ii]] = ii;
  }
  RNG rng_to_generate_swaps;
  auto nodes_copy = nodes;
  rng_to_generate_swaps.do_shuffle(nodes_copy);
  const auto node_final_positions = nodes_copy;

  problem_ss << " Node mapping:";
  NodeMapping node_mapping;
  for (size_t ii = 0; ii < nodes.size(); ++ii) {
    problem_ss << "\ni=" << ii << " : " << node_final_positions[ii].repr()
               << " -> " << nodes[ii].repr();
    node_mapping[node_final_positions[ii]] = nodes[ii];
  }
  CHECK(
      problem_ss.str() ==
      "24 nodes; 46 edges. Node mapping:\n"
      "i=0 : gridNode[0, 0, 0] -> gridNode[0, 0, 0]\n"
      "i=1 : gridNode[0, 3, 0] -> gridNode[0, 0, 1]\n"
      "i=2 : gridNode[2, 1, 0] -> gridNode[0, 1, 0]\n"
      "i=3 : gridNode[0, 1, 1] -> gridNode[0, 1, 1]\n"
      "i=4 : gridNode[2, 2, 0] -> gridNode[0, 2, 0]\n"
      "i=5 : gridNode[1, 1, 1] -> gridNode[0, 2, 1]\n"
      "i=6 : gridNode[0, 0, 1] -> gridNode[0, 3, 0]\n"
      "i=7 : gridNode[0, 3, 1] -> gridNode[0, 3, 1]\n"
      "i=8 : gridNode[1, 3, 0] -> gridNode[1, 0, 0]\n"
      "i=9 : gridNode[1, 0, 0] -> gridNode[1, 0, 1]\n"
      "i=10 : gridNode[2, 2, 1] -> gridNode[1, 1, 0]\n"
      "i=11 : gridNode[0, 1, 0] -> gridNode[1, 1, 1]\n"
      "i=12 : gridNode[2, 0, 1] -> gridNode[1, 2, 0]\n"
      "i=13 : gridNode[1, 2, 1] -> gridNode[1, 2, 1]\n"
      "i=14 : gridNode[1, 3, 1] -> gridNode[1, 3, 0]\n"
      "i=15 : gridNode[1, 0, 1] -> gridNode[1, 3, 1]\n"
      "i=16 : gridNode[2, 0, 0] -> gridNode[2, 0, 0]\n"
      "i=17 : gridNode[2, 1, 1] -> gridNode[2, 0, 1]\n"
      "i=18 : gridNode[0, 2, 1] -> gridNode[2, 1, 0]\n"
      "i=19 : gridNode[1, 2, 0] -> gridNode[2, 1, 1]\n"
      "i=20 : gridNode[0, 2, 0] -> gridNode[2, 2, 0]\n"
      "i=21 : gridNode[1, 1, 0] -> gridNode[2, 2, 1]\n"
      "i=22 : gridNode[2, 3, 0] -> gridNode[2, 3, 0]\n"
      "i=23 : gridNode[2, 3, 1] -> gridNode[2, 3, 1]");

  // Calculate swaps to enact the permutation.
  const auto node_swaps = get_swaps(arch, node_mapping);
  // This will hopefully decrease over time
  // as we improve the algorithm.
  CHECK(node_swaps.size() == 29);

  // Go back to the original configuration, and perform the swaps.
  nodes_copy = nodes;
  for (const auto& node_swap : node_swaps) {
    REQUIRE(allowed_edges_map.at(node_swap.first).count(node_swap.second) != 0);
    const auto index1 = original_vertex_indices.at(node_swap.first);
    const auto index2 = original_vertex_indices.at(node_swap.second);
    REQUIRE(index1 != index2);
    std::swap(nodes_copy[index1], nodes_copy[index2]);
  }
  REQUIRE(nodes_copy == node_final_positions);
}

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket

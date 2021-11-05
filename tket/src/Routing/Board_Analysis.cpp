// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "Routing.hpp"

namespace tket {

bool node_active(const qubit_bimap_t& map, Node node) {
  const bool found = map.right.find(node) != map.right.end();
  return found;
}

Node Routing::find_best_inactive_node(
    const Node& target_node, const Architecture& arc) const {
  const unsigned diameter = arc.get_diameter();
  for (unsigned k = 1; k <= diameter; k++) {
    std::vector<Node> potential_nodes = arc.nodes_at_distance(target_node, k);
    for (Node potential : potential_nodes) {
      if (!node_active(qmap, potential)) {
        return potential;
      }
    }
  }
  throw ArchitectureFull();  // gotta hope you never get here...
}

void Routing::activate_node(const Node& node) {
  current_arc_.add_node(node);
  for (Node neigh : original_arc_.get_neighbour_nodes(node)) {
    if (node_active(qmap, neigh)) {
      if (original_arc_.edge_exists(node, neigh)) {
        current_arc_.add_connection(node, neigh);
      }
      if (original_arc_.edge_exists(neigh, node)) {
        current_arc_.add_connection(neigh, node);
      }
    }
  }
}

void Routing::reactivate_qubit(const Qubit& qb, const Qubit& target) {
  // finds 'best' available node
  Node node = find_best_inactive_node(qmap.left.at(target), original_arc_);

  // updates qmap and initial maps to reflect this qb being at that node
  activate_node(node);
  std::pair<Qubit, Node> new_in = {qb, node};
  qmap.left.insert(new_in);
  init_map.left.insert(new_in);
}

}  // namespace tket

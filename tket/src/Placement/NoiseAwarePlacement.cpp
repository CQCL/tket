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

#include "Placement/Placement.hpp"
#include "Utils/HelperFunctions.hpp"

namespace tket {

double NoiseAwarePlacement::cost_placement(
    const boost::bimap<Qubit, Node>& map, const Circuit& circ_,
    const QubitGraph& q_graph) const {
  double cost = 0.0;
  if (circ_.n_gates() == 0 || circ_.n_qubits() == 0) {
    return cost;
  }
  const int approx_depth = circ_.n_gates() / circ_.n_qubits() + 1;
  // constants for scaling single qubit error
  constexpr double c1 = 0.5;
  constexpr double d1 = 1 - 1 / c1;
  for (auto [qb, node] : map) {
    node_set_t neighbours = this->architecture_.get_neighbour_nodes(node);
    double edge_sum = 1.0;
    for (const Node& neighbour : neighbours) {
      double fwd_edge_weighting = 1.0, bck_edge_weighting = 1.0;
      // check if neighbour node is mapped
      auto nei_qb_it = map.right.find(neighbour);
      if (nei_qb_it == map.right.end()) continue;

      auto place_interactions_boost = [&](unsigned edge_v) {
        return -1 * (edge_v - this->maximum_pattern_depth_ );
      };
      
      // check if either directed interaction exists
      // if edge is used by interaction in mapping, weight edge higher

      unsigned edge_val = q_graph.get_connection_weight(qb, nei_qb_it->second);
      if (edge_val) {
        fwd_edge_weighting += place_interactions_boost(edge_val);
      } else {
        // edge_exists = q_graph.edge_exists(nei_qb_it->second, qb);
        edge_val = q_graph.get_connection_weight(nei_qb_it->second, qb);
        if (edge_val) {
          bck_edge_weighting += place_interactions_boost(edge_val);
        }
      }
      gate_error_t fwd_error =
          this->characterisation_.get_error({node, neighbour});
      gate_error_t bck_error =
          this->characterisation_.get_error({neighbour, node});
      if(fwd_error){
        edge_sum += fwd_edge_weighting * fwd_error;
      }
      if(bck_error){
        edge_sum += bck_edge_weighting * bck_error;
      }
    }
    // bigger edge sum -> smaller cost
    cost += 1.0 / (edge_sum);
    // add error rate of node
    gate_error_t single_error = this->characterisation_.get_error(node);
    if(single_error){
      cost += d1 + 1.0 / ((1.0 - single_error) + c1);
    }
    readout_error_t readout_error =
        this->characterisation_.get_readout_error(node);
    if(readout_error){
      cost += (d1 + 1.0 / ((1.0 - readout_error) + c1)) / (approx_depth * 20);
    }
  }
  return cost;
}

std::vector<std::map<Qubit, Node>> NoiseAwarePlacement::rank_maps(
    const std::vector<boost::bimap<Qubit, Node>>& placement_maps,
    const Circuit& circ_,
    const std::vector<GraphPlacement::WeightedEdge>& pattern_edges) const {
  // each map is costed, sorted and returned
  std::vector<std::map<Qubit, Node>> return_placement_maps;
  std::vector<double> costs;
  QubitGraph q_graph =
      this->construct_pattern_graph(pattern_edges, circ_.n_qubits());
  for (const boost::bimap<Qubit, Node>& map : placement_maps) {
    double cost = this->cost_placement(map, circ_, q_graph);
    unsigned size = costs.size();
    for (unsigned i = 0; i < costs.size(); i++) {
      if (cost > costs[i]) {
        costs.insert(costs.begin() + i, cost);
      std::map<Qubit, Node> map_qn = bimap_to_map(map.left);
        return_placement_maps.insert(
            return_placement_maps.begin() + i, bimap_to_map(map.left));
        break;
      }
    }
    if (size == costs.size()) {
      costs.push_back(cost);
      std::map<Qubit, Node> map_qn = bimap_to_map(map.left);
      return_placement_maps.push_back(bimap_to_map(map.left));
    }
  }
  return return_placement_maps;
}

std::vector<std::map<Qubit, Node>> NoiseAwarePlacement::get_all_placement_maps(
    const Circuit& circ_, unsigned matches) {
  std::vector<WeightedEdge> weighted_pattern_edges =
      this->default_pattern_weighting(circ_);
  std::vector<boost::bimap<Qubit, Node>> placement_maps =
      this->get_all_weighted_subgraph_monomorphisms(
          circ_, weighted_pattern_edges, true);
  std::cout << "n placement maps: " << placement_maps.size() << std::endl;
  std::vector<std::map<Qubit, Node>> ranked_placement_maps =
      this->rank_maps(placement_maps, circ_, weighted_pattern_edges);
  if (ranked_placement_maps.size() > matches) {
    ranked_placement_maps.resize(matches);
  }
  return ranked_placement_maps;
}

}  // namespace tket
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

#include "tket/Placement/Placement.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {

NoiseAwarePlacement::NoiseAwarePlacement(
    const Architecture& _architecture,
    std::optional<avg_node_errors_t> _node_errors,
    std::optional<avg_link_errors_t> _link_errors,
    std::optional<avg_readout_errors_t> _readout_errors,
    unsigned _maximum_matches, unsigned _timeout,
    unsigned _maximum_pattern_gates, unsigned _maximum_pattern_depth)
    : GraphPlacement(
          _architecture, _maximum_matches, _timeout, _maximum_pattern_gates,
          _maximum_pattern_depth) {
  architecture_ = _architecture;
  this->weighted_target_edges = this->default_target_weighting(architecture_);
  this->extended_target_graphs = {
      this->construct_target_graph(weighted_target_edges, 0)
          .get_undirected_connectivity()};
  characterisation_ = {
      _node_errors ? *_node_errors : avg_node_errors_t(),
      _link_errors ? *_link_errors : avg_link_errors_t(),
      _readout_errors ? *_readout_errors : avg_readout_errors_t()};
}

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
        return this->maximum_pattern_depth_ - edge_v + 1;
      };

      // check if either directed interaction exists
      // if edge is used by interaction in mapping, weight edge higher
      unsigned edge_val = q_graph.get_connection_weight(qb, nei_qb_it->second);
      if (edge_val) {
        fwd_edge_weighting += place_interactions_boost(edge_val);
      } else {
        edge_val = q_graph.get_connection_weight(nei_qb_it->second, qb);
        if (edge_val) {
          bck_edge_weighting += place_interactions_boost(edge_val);
        }
      }
      gate_error_t fwd_error =
          this->characterisation_.get_error({node, neighbour});
      gate_error_t bck_error =
          this->characterisation_.get_error({neighbour, node});

      if (fwd_error < 1.0 && bck_error < 1.0) {
        edge_sum += fwd_edge_weighting * (1.0 - fwd_error);
        edge_sum += bck_edge_weighting * (1.0 - bck_error);
      }
    }
    // bigger edge sum -> smaller cost
    cost += 1.0 / (edge_sum);
    // add error rate of node
    gate_error_t single_error = this->characterisation_.get_error(node);
    cost += d1 + 1.0 / ((1.0 - single_error) + c1);
    readout_error_t readout_error =
        this->characterisation_.get_readout_error(node);
    if (readout_error) {
      cost += (d1 + 1.0 / ((1.0 - readout_error) + c1)) / (approx_depth * 20);
    }
  }
  return cost;
}

std::vector<boost::bimap<Qubit, Node>> NoiseAwarePlacement::rank_maps(
    const std::vector<boost::bimap<Qubit, Node>>& placement_maps,
    const Circuit& circ_,
    const std::vector<GraphPlacement::WeightedEdge>& pattern_edges) const {
  // each map is costed, sorted and returned
  // only equal best costed maps are returned
  std::vector<boost::bimap<Qubit, Node>> return_placement_maps;
  double best_cost = 0;
  QubitGraph q_graph =
      this->construct_pattern_graph(pattern_edges, circ_.n_qubits());
  for (const boost::bimap<Qubit, Node>& map : placement_maps) {
    double cost = this->cost_placement(map, circ_, q_graph);
    if (return_placement_maps.empty() || cost < best_cost) {
      best_cost = cost;
      return_placement_maps = {map};
    } else if (cost == best_cost) {
      return_placement_maps.push_back(map);
    }
  }
  return return_placement_maps;
}

std::vector<std::map<Qubit, Node>> NoiseAwarePlacement::get_all_placement_maps(
    const Circuit& circ_, unsigned matches) const {
  std::vector<WeightedEdge> weighted_pattern_edges =
      this->default_pattern_weighting(circ_);

  if (weighted_pattern_edges.empty()) {
    // => there are no two-qubit gates in the input circuit
    // in this case, as this method is "noise-aware", assign
    // qubits in the circuit to the lowest-error architecture nodes
    std::vector<std::pair<gate_error_t, Node>> all_node_errors;
    for (const Node& node : this->architecture_.nodes()) {
      all_node_errors.push_back(
          {this->characterisation_.get_error(node), node});
    }
    std::sort(
        all_node_errors.begin(), all_node_errors.end(),
        [](const auto& lhs, const auto& rhs) {
          // Compare based on the float values in descending order
          return lhs.first < rhs.first;
        });
    std::vector<Qubit> circ_qubits = circ_.all_qubits();
    TKET_ASSERT(all_node_errors.size() >= circ_qubits.size());
    std::map<Qubit, Node> placement_map;
    for (std::size_t i = 0; i < circ_qubits.size(); i++) {
      placement_map.insert({circ_qubits[i], all_node_errors[i].second});
    }
    return {placement_map};
  }

  std::vector<boost::bimap<Qubit, Node>> placement_maps =
      this->get_all_weighted_subgraph_monomorphisms(
          circ_, weighted_pattern_edges, true);
  std::vector<boost::bimap<Qubit, Node>> ranked_placement_maps =
      this->rank_maps(placement_maps, circ_, weighted_pattern_edges);
  std::vector<std::map<Qubit, Node>> all_qmaps;
  QubitGraph::UndirectedConnGraph pattern_graph =
      this->construct_pattern_graph(weighted_pattern_edges, circ_.n_qubits())
          .get_undirected_connectivity();
  for (unsigned i = 0; i < ranked_placement_maps.size() && i < matches; i++) {
    all_qmaps.push_back(convert_bimap(ranked_placement_maps[i], pattern_graph));
  }
  return all_qmaps;
}

DeviceCharacterisation NoiseAwarePlacement::get_characterisation() const {
  return this->characterisation_;
}

void NoiseAwarePlacement::set_characterisation(
    const DeviceCharacterisation& characterisation) {
  this->characterisation_ = characterisation;
}

}  // namespace tket
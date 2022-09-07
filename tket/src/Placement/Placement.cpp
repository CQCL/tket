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

bool Placement::place(Circuit& circ_) const {
  std::map<UnitID, UnitID> map_ = this->get_placement_map(circ_);
  std::map<Qubit, Node> recast_map;
  for (const auto& entry : map_) {
    recast_map.insert({Qubit(entry.first), Node(entry.second)});
  }
  return circ_.rename_units(recast_map);
}

std::map<UnitID, UnitID> Placement::get_placement_map(
    const Circuit& circ_) const {
  std::vector<std::map<UnitID, UnitID>> all_maps =
      this->get_all_placement_maps(circ_);
  // basic handling to avoid segmentation faults, as placement method may not
  // return any valid map
  auto it = all_maps.begin();
  if (it != all_maps.end()) {
    return *it;
  } else {
    return {};
  }
}

std::vector<std::map<UnitID, UnitID>> Placement::get_all_placement_maps(
    const Circuit& /*circ_*/) const {
  return {};
}

const std::vector<GraphPlacement::WeightedEdge>
GraphPlacement::default_weighting(const Circuit& circuit) {
  Placement::Frontier frontier(circuit);
  unsigned max_gates = 100, max_depth = 100, gate_counter = 0;
  std::vector<GraphPlacement::WeightedEdge> weights;
  for (unsigned i = 0;
       i < max_depth && gate_counter < max_gates && !frontier.slice->empty();
       i++) {
    for (const Vertex& vert : *frontier.slice) {
      EdgeVec q_out_edges =
          circuit.get_out_edges_of_type(vert, EdgeType::Quantum);
      unsigned n_q_edges = q_out_edges.size();
      if (n_q_edges == 2) {
        UnitID uid_0, uid_1;
        Edge edge_0 = q_out_edges[0];
        Edge edge_1 = q_out_edges[1];
        bool match_0 = false, match_1 = false;
        for (const std::pair<UnitID, Edge>& pair :
             frontier.quantum_out_edges->get<TagKey>()) {
          if (!match_0 && pair.second == edge_0) {
            uid_0 = pair.first;
          }
          if (!match_1 && pair.second == edge_1) {
            uid_1 = pair.first;
          }
          if (match_0 && match_1) {
            break;
          }
        }

        bool match_weight = false;
        for (WeightedEdge& weighted_edge : weights) {
          //   std::pair<UnitID, UnitID> uids = uid_w.first;
          if ((weighted_edge.node0 == uid_0 && weighted_edge.node1 == uid_1) ||
              (weighted_edge.node0 == uid_1 && weighted_edge.node1 == uid_0)) {
            // actually update the weight here, i.e. this is the "magic"
            match_weight = true;
            weighted_edge.weight += (max_depth - i);
            break;
          }
        }
        if (!match_weight) {
          weights.push_back({uid_0, uid_1, max_depth - i});
        }
        gate_counter++;
      }
      if (n_q_edges > 2) {
        throw std::invalid_argument(
            "Can only weight for Circuits with maximum two qubit gates.");
      }
    }
  }
  return weights;
}

QubitGraph GraphPlacement::construct_pattern_graph(
    const std::vector<WeightedEdge>& edges) const {
  QubitGraph q_graph;
  for (const WeightedEdge& weighted_edge : edges) {
    Node node0 = Node(weighted_edge.node0);
    Node node1 = Node(weighted_edge.node1);
    bool e_01_exists = q_graph.edge_exists(node0, node1);
    bool e_10_exists = q_graph.edge_exists(node1, node0);
    if (e_01_exists || e_10_exists) {
      throw std::invalid_argument(
          "Graph can only have on edge between a pair of Node.");
    }
    if (!q_graph.node_exists(node0)) {
      q_graph.add_node(node0);
    }
    if (!q_graph.node_exists(node1)) {
      q_graph.add_node(node1);
    }
    if (weighted_edge.weight > 0) {
      q_graph.add_connection(node0, node1, weighted_edge.weight);
    }
  }
  return q_graph;
}

std::vector<std::map<UnitID, UnitID>> GraphPlacement::get_all_placement_maps(
    const Circuit& circ_) const {
  std::vector<WeightedEdge> weighted_edges = this->weight_pattern_graph_(circ_);
  QubitGraph pattern_graph = this->construct_pattern_graph(weighted_edges);
  Architecture::UndirectedConnGraph target_graph =
      this->arc_.get_undirected_connectivity();
  std::vector<boost::bimap<Qubit, Node>> all_bimaps =
      this->get_weighted_subgraph_monomorphisms(
          pattern_graph, target_graph, this->maximum_matches_, this->timeout_);
  std::vector<std::map<Qubit, Node>> all_qmaps;
  for (boost::bimap<Qubit, Node>& bm : all_bimaps) {
    all_qmaps.push_back(bimap_to_map(bm.left));
  }
  //   TODO: Update this to convert to UnitID to start with, this is
  //   unncessary...
  std::vector<std::map<UnitID, UnitID>> all_uidmap;
  for (const auto& map : all_qmaps) {
    std::map<UnitID, UnitID> uid_map;
    for (const auto& q_n : map) {
      uid_map.insert({UnitID(q_n.first), UnitID(q_n.second)});
    }
  }
  return all_uidmap;
}
}  // namespace tket
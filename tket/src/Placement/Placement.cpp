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

bool Placement::place(
    Circuit& circ_, std::shared_ptr<unit_bimaps_t> compilation_map) {
  if (circ_.n_qubits() > this->architecture_.n_nodes()) {
    throw std::invalid_argument(
        "Circuit has more qubits than Architecture has nodes.");
  }
  std::map<Qubit, Node> map_ = this->get_placement_map(circ_);
  return this->place_with_map(circ_, map_, compilation_map);
}

bool Placement::place_with_map(
    Circuit& circ_, std::map<Qubit, Node>& map_,
    std::shared_ptr<unit_bimaps_t> compilation_map) {
  bool changed = circ_.rename_units(map_);
  changed |= update_maps(compilation_map, map_, map_);
  return changed;
}

std::map<Qubit, Node> Placement::get_placement_map(const Circuit& circ_) {
  std::vector<std::map<Qubit, Node>> all_maps =
      this->get_all_placement_maps(circ_, 1);
  // basic handling to avoid segmentation faults, as placement method may not
  // return any valid map
  auto it = all_maps.begin();
  if (it != all_maps.end()) {
    return *it;
  } else {
    return {};
  }
}

std::vector<std::map<Qubit, Node>> Placement::get_all_placement_maps(
    const Circuit& circ_, unsigned /*matches*/) {
  std::map<Qubit, Node> placement;
  qubit_vector_t to_place;
  std::vector<Node> placed;

  // Find which/if any qubits need placing
  for (const Qubit& q : circ_.all_qubits()) {
    Node n(q);
    if (!this->architecture_.node_exists(n)) {
      to_place.push_back(n);
    } else {
      placed.push_back(n);
      // if already placed, make sure qubit retains placement
      placement.insert({n, n});
    }
  }
  // avoid doing std::set_difference unless qubits need to be placed
  unsigned n_placed = to_place.size();
  if (n_placed > 0) {
    std::vector<Node> difference,
        architecture_nodes = this->architecture_.get_all_nodes_vec();
    std::set_difference(
        architecture_nodes.begin(), architecture_nodes.end(), placed.begin(),
        placed.end(), std::inserter(difference, difference.begin()));
    // should always be enough remaining qubits to assign unplaced qubits to
    if (difference.size() < n_placed) {
      throw std::invalid_argument(
          "There are more unplaced Qubit in Circuit then there are free Nodes "
          "in Architecture.");
    }
    for (unsigned i = 0; i < n_placed; i++) {
      // naively assign each qubit to some free node
      placement.insert({to_place[i], difference[i]});
    }
  }
  return {placement};
}

const std::vector<GraphPlacement::WeightedEdge>
GraphPlacement::default_pattern_weighting(const Circuit& circuit) {
  GraphPlacement::Frontier frontier(circuit);
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
          if ((weighted_edge.node0 == uid_0 && weighted_edge.node1 == uid_1) ||
              (weighted_edge.node0 == uid_1 && weighted_edge.node1 == uid_0)) {
            // actually update the weight here, i.e. this is the "magic"
            match_weight = true;
            weighted_edge.weight += unsigned(max_depth - i);
            break;
          }
        }
        if (!match_weight) {
          weights.push_back({uid_0, uid_1, unsigned(max_depth - i), 0});
        }
        gate_counter++;
      }
      if (n_q_edges > 2) {
        throw std::invalid_argument(
            "Can only weight for Circuits with maximum two qubit gates.");
      }
    }
    frontier.next_slicefrontier();
  }
  return weights;
}

const std::vector<GraphPlacement::WeightedEdge>
GraphPlacement::default_target_weighting(Architecture& passed_architecture) {
  std::vector<Node> all_nodes = passed_architecture.get_all_nodes_vec();
  std::vector<GraphPlacement::WeightedEdge> weights;
  unsigned diameter = passed_architecture.get_diameter();
  auto it = all_nodes.begin();
  while (it != all_nodes.end()) {
    auto jt = it;
    ++jt;
    while (jt != all_nodes.end()) {
      unsigned distance = passed_architecture.get_distance(*it, *jt);
      weights.push_back(
          {*it, *jt, unsigned(diameter + 1 - distance), distance});
      ++jt;
    }
    ++it;
  }
  return weights;
}

QubitGraph GraphPlacement::construct_pattern_graph(
    const std::vector<WeightedEdge>& edges, unsigned max_out_degree) const {
  QubitGraph q_graph;

  for (const WeightedEdge& weighted_edge : edges) {
    Node node0 = Node(weighted_edge.node0);
    Node node1 = Node(weighted_edge.node1);
    if (!q_graph.node_exists(node0)) {
      q_graph.add_node(node0);
    }
    if (!q_graph.node_exists(node1)) {
      q_graph.add_node(node1);
    }
    bool e_01_exists = q_graph.edge_exists(node0, node1);
    bool e_10_exists = q_graph.edge_exists(node1, node0);
    if (e_01_exists || e_10_exists) {
      throw std::invalid_argument(
          "Graph can only have a single edge between a pair of Node.");
    }
    if (weighted_edge.weight > 0 &&
        q_graph.get_out_degree(node0) < max_out_degree &&
        q_graph.get_out_degree(node1) < max_out_degree) {
      q_graph.add_connection(node0, node1, weighted_edge.weight);
    }
  }
  return q_graph;
}

Architecture GraphPlacement::construct_target_graph(
    const std::vector<WeightedEdge>& edges, unsigned distance) const {
  Architecture architecture;
  for (const WeightedEdge& weighted_edge : edges) {
    Node node0 = Node(weighted_edge.node0);
    Node node1 = Node(weighted_edge.node1);
    if (!architecture.node_exists(node0)) {
      architecture.add_node(node0);
    }
    if (!architecture.node_exists(node1)) {
      architecture.add_node(node1);
    }
    bool e_01_exists = architecture.edge_exists(node0, node1);
    bool e_10_exists = architecture.edge_exists(node1, node0);
    if (e_01_exists || e_10_exists) {
      throw std::invalid_argument(
          "Graph can only have a single edge between a pair of Node.");
    }
    if (weighted_edge.weight > 0 && weighted_edge.distance <= distance + 1) {
      architecture.add_connection(node0, node1, weighted_edge.weight);
    }
  }
  return architecture;
}

std::vector<std::map<Qubit, Node>> GraphPlacement::get_all_placement_maps(
    const Circuit& circ_, unsigned matches) {
  if (circ_.n_qubits() > this->architecture_.n_nodes()) {
    throw std::invalid_argument(
        "Circuit has more qubits than Architecture has nodes.");
  }
  unsigned n_qubits = circ_.n_qubits();
  if (n_qubits == 0) {
    return {{}};
  }

  /** The weighted subgraph monomorphism tool from TK-WSM is efficient at
   * returning nothing when no subgraph monomorphism can be found. The otherside
   * to this is that it typically finds no "partial" solutions.
   *
   * Therefore, to provide "good" program to physical qubit assignments we
   * must emulate finding partial solutions with the TK-WSM.
   *
   * At GraphPlacement object construction potential target graph edges are
   * weighted from the given Architecture At get_all_placement_maps calls
   * potential pattern graph edges are weighted from the given Circuit
   *
   * From these edges, weighted pattern graphs and weighted target graphs are
   * constructed until solutions are found. The approach is to move from optimal
   * solutions to partial assumptive solutions.
   *
   * As we know TK-WSM will quickly return false if a subgraph monomorphism is
   * impossible to find, so we can use it to build pattern and target graphs
   * that are valid.
   *
   * As finding the distance between all pairs of Nodes in an Architecture is
   * expensive, we cache the constructed target graphs.
   *
   *
   *
   * Finally, for each returned assignment, we remove assignments that do not
   * have interactions with all their neighbours. This allows Routing to
   * dynamically assign them at point they are encountered.
   *
   * Also note that given the symmetry of typical architecture graphs, at the
   * point a solution is found there are often many valid assignments.
   */

  std::vector<WeightedEdge> weighted_pattern_edges =
      this->weight_pattern_graph_(circ_);

  if (weighted_pattern_edges.empty()) {
    return {{}};
  }

  // We store pattern graphs as they're constructed, and check each of them in
  // less complex order when a new target graph is constructed
  std::vector<QubitGraph::UndirectedConnGraph> all_pattern_graphs;
  std::vector<boost::bimap<Qubit, Node>> all_bimaps;
  unsigned incrementer = 0;
  while (all_bimaps.empty()) {
    /**
     * Note that this is the while loop condition as this will always terminate
     * As eventually an edge will be added between every Node on the
     * Architecture, meaning a solution will be found.
     */
    if (extended_target_graphs.size() <= incrementer) {
      extended_target_graphs.push_back(
          this->construct_target_graph(this->weighted_target_edges, incrementer)
              .get_undirected_connectivity());
      TKET_ASSERT(extended_target_graphs.size() - 1 == incrementer);
    }
    TKET_ASSERT(extended_target_graphs.size() > incrementer);

    // For each increment we construct a smaller pattern graph
    all_pattern_graphs.push_back(
        this->construct_pattern_graph(
                weighted_pattern_edges, n_qubits - incrementer - 1)
            .get_undirected_connectivity());
    // For each pattern graph constructed, we attempt to find
    // a subgraph monomorphism for the new target graph
    // From more full to elss full
    auto it = all_pattern_graphs.begin();
    while (it != all_pattern_graphs.end() && all_bimaps.empty()) {
      all_bimaps = get_weighted_subgraph_monomorphisms(
          *it, extended_target_graphs[incrementer], this->maximum_matches_,
          this->timeout_);
      ++it;
    }
    incrementer++;
  }

  std::vector<std::map<Qubit, Node>> all_qmaps;
  unsigned counter = 0;
  for (auto it = all_bimaps.begin();
       it != all_bimaps.end() && counter < matches; ++it) {
    // TODO: clean up the solution by removing low cost nodes/unconnected nodes
    all_qmaps.push_back(bimap_to_map(it->left));
    ++counter;
  }
  return all_qmaps;
}

void to_json(nlohmann::json& j, const Placement::Ptr& placement_ptr) {
  j["architecture"] = placement_ptr->get_architecture_ref();
  if (std::shared_ptr<GraphPlacement> cast_placer =
          std::dynamic_pointer_cast<GraphPlacement>(placement_ptr)) {
    j["type"] = "GraphPlacement";
    j["matches"] = cast_placer->get_maximum_matches();
    j["timeout"] = cast_placer->get_timeout();
  } else {
    j["type"] = "Placement";
  }
}

void from_json(const nlohmann::json& j, Placement::Ptr& placement_ptr) {
  std::string classname = j.at("type").get<std::string>();
  Architecture arc = j.at("architecture").get<Architecture>();
  if (classname == "GraphPlacement") {
    unsigned matches = j.at("matches").get<unsigned>();
    unsigned timeout = j.at("timeout").get<unsigned>();
    placement_ptr = std::make_shared<GraphPlacement>(arc, matches, timeout);
  } else {
    placement_ptr = std::make_shared<Placement>(arc);
  }
}
}  // namespace tket
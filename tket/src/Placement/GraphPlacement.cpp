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

#include <chrono>

#include "tket/Placement/Placement.hpp"

typedef std::chrono::steady_clock Clock;

namespace tket {

GraphPlacement::GraphPlacement(
    const Architecture& _architecture, unsigned _maximum_matches,
    unsigned _timeout, unsigned _maximum_pattern_gates,
    unsigned _maximum_pattern_depth)

    : maximum_matches_(_maximum_matches),
      timeout_(_timeout),
      maximum_pattern_gates_(_maximum_pattern_gates),
      maximum_pattern_depth_(_maximum_pattern_depth) {
  architecture_ = _architecture;
  this->weighted_target_edges = this->default_target_weighting(architecture_);
  this->extended_target_graphs = {
      this->construct_target_graph(weighted_target_edges, 0)
          .get_undirected_connectivity()};
}

const std::vector<GraphPlacement::WeightedEdge>
GraphPlacement::default_pattern_weighting(const Circuit& circuit) const {
  GraphPlacement::Frontier frontier(circuit);
  unsigned gate_counter = 0;
  std::vector<GraphPlacement::WeightedEdge> weights;
  for (unsigned i = 0;
       i < this->maximum_pattern_depth_ &&
       gate_counter < this->maximum_pattern_gates_ && !frontier.slice->empty();
       i++) {
    for (const Vertex& vert : *frontier.slice) {
      if (circuit.get_OpType_from_Vertex(vert) == OpType::Barrier) continue;
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
            weighted_edge.weight += unsigned(this->maximum_pattern_depth_ - i);
            break;
          }
        }
        if (!match_weight) {
          weights.push_back(
              {uid_0, uid_1, unsigned(this->maximum_pattern_depth_ - i), 0});
        }
        gate_counter++;
      }
      if (n_q_edges > 2) {
        throw std::invalid_argument(
            "Can only weight for Circuits with maximum two qubit quantum "
            "gates.");
      }
    }
    frontier.next_slicefrontier();
  }
  return weights;
}

const std::vector<GraphPlacement::WeightedEdge>
GraphPlacement::default_target_weighting(
    Architecture& passed_architecture) const {
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
  // Note that edges are in weight order, so this adds higher weight edges first
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

std::vector<boost::bimap<Qubit, Node>>
GraphPlacement::get_all_weighted_subgraph_monomorphisms(
    const Circuit& circ_,
    const std::vector<GraphPlacement::WeightedEdge>& weighted_pattern_edges,
    bool return_best) const {
  // we routinely check time, and throw a runtime_error if breached
  const auto init_start = Clock::now();
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

  if (weighted_pattern_edges.empty()) {
    return {{}};
  }

  if (std::chrono::duration_cast<std::chrono::milliseconds>(
          Clock::now() - init_start)
          .count() > this->timeout_) {
    throw std::runtime_error(
        "GraphPlacement execution time has exceeded allowed limits.");
  }
  // We store pattern graphs as they're constructed, and check each of them in
  // less complex order when a new target graph is constructed
  std::vector<QubitGraph::UndirectedConnGraph> all_pattern_graphs;
  std::vector<boost::bimap<Qubit, Node>> all_bimaps;
  unsigned incrementer = 0, last_edges = 0;
  while (all_bimaps.empty()) {
    // we check timeout not reached regularly
    if (std::chrono::duration_cast<std::chrono::milliseconds>(
            Clock::now() - init_start)
            .count() > this->timeout_) {
      throw std::runtime_error(
          "GraphPlacement execution time has exceeded allowed limits.");
    }
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
    QubitGraph::UndirectedConnGraph pattern_graph =
        this->construct_pattern_graph(
                weighted_pattern_edges, n_qubits - incrementer - 1)
            .get_undirected_connectivity();
    // It's possible that no edges are removed, so only add new graph if it
    // has a different number of edges (i.e. is different)
    unsigned n_edges = boost::num_edges(pattern_graph);
    if (last_edges != n_edges) {
      all_pattern_graphs.push_back(pattern_graph);
      last_edges = n_edges;
    }
    // For each pattern graph constructed, we attempt to find
    // a subgraph monomorphism for the new target graph
    // From more full to less full
    auto it = all_pattern_graphs.begin();
    while (it != all_pattern_graphs.end() && all_bimaps.empty()) {
      all_bimaps = get_weighted_subgraph_monomorphisms(
          *it, extended_target_graphs[incrementer], this->maximum_matches_,
          this->timeout_ -
              std::chrono::duration_cast<std::chrono::milliseconds>(
                  Clock::now() - init_start)
                  .count(),
          return_best);

      ++it;
    }
    incrementer++;
  }
  return all_bimaps;
}

std::map<Qubit, Node> GraphPlacement::convert_bimap(
    boost::bimap<Qubit, Node>& bimap,
    const QubitGraph::UndirectedConnGraph& pattern_graph) const {
  /**
   * For each assignment, find Qubit on adjacent Node
   * If a majority of adjacent Qubits don't interacting with
   * assignment, then remove assignment
   */
  std::map<Qubit, Node> out_map;
  // construct a map between pattern graph Qubit and Vertex
  std::map<Qubit, unsigned> qubit_vertex_map;
  auto vertex_iter_pair = boost::vertices(pattern_graph);
  while (vertex_iter_pair.first != vertex_iter_pair.second) {
    unsigned vert = *vertex_iter_pair.first;
    qubit_vertex_map.insert({pattern_graph[vert], vert});
    ++vertex_iter_pair.first;
  }
  for (const auto& entry : bimap.left) {
    auto qvm_it = qubit_vertex_map.find(entry.first);
    TKET_ASSERT(qvm_it != qubit_vertex_map.end());
    unsigned entry_vertex = qvm_it->second;
    unsigned n_pattern_edges = boost::out_degree(entry_vertex, pattern_graph);
    std::set<Node> neighbour_nodes =
        this->architecture_.get_neighbour_nodes(entry.second);
    unsigned n_target_edges = neighbour_nodes.size();
    unsigned n_interacting = 0;
    for (const Node& node : neighbour_nodes) {
      auto it = bimap.right.find(node);
      // Node may be not be assigned to
      // If it is, check if Qubit are interacting
      // If not, decrement n_interacting
      if (it != bimap.right.end()) {
        auto qvm_it = qubit_vertex_map.find(it->second);
        TKET_ASSERT(qvm_it != qubit_vertex_map.end());
        auto [_, exists] =
            boost::edge(entry_vertex, qvm_it->second, pattern_graph);
        if (exists) {
          n_interacting++;
        }
      }
    }
    if (std::min(n_pattern_edges, n_target_edges) - n_interacting <=
        n_interacting) {
      out_map.insert({entry.first, entry.second});
    }
  }
  return out_map;
}

std::vector<std::map<Qubit, Node>> GraphPlacement::get_all_placement_maps(
    const Circuit& circ_, unsigned matches) const {
  std::vector<WeightedEdge> weighted_pattern_edges =
      this->default_pattern_weighting(circ_);
  std::vector<boost::bimap<Qubit, Node>> all_bimaps =
      this->get_all_weighted_subgraph_monomorphisms(
          circ_, weighted_pattern_edges, false);
  std::vector<std::map<Qubit, Node>> all_qmaps;
  QubitGraph::UndirectedConnGraph pattern_graph =
      this->construct_pattern_graph(weighted_pattern_edges, circ_.n_qubits())
          .get_undirected_connectivity();
  for (unsigned i = 0; i < all_bimaps.size() && i < matches; i++) {
    all_qmaps.push_back(convert_bimap(all_bimaps[i], pattern_graph));
  }
  return all_qmaps;
}
}  // namespace tket
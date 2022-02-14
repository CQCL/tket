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

//#define DEBUG
#include <algorithm>
#include <numeric>
#include <optional>
#include <queue>

#include "Architecture/Architecture.hpp"
#include "Graphs/Utils.hpp"
#include "Placement.hpp"

namespace tket {

std::set<Qubit> interacting_qbs(const Circuit& circ) {
  std::set<Qubit> qbs;
  for (const Qubit& qb : circ.all_qubits()) {
    Edge e = circ.get_nth_out_edge(circ.get_in(qb), 0);
    Vertex terminal = circ.target(circ.skip_irrelevant_edges(e));
    if (!circ.detect_final_Op(terminal)) {
      qbs.insert(qb);
    }
  }

  return qbs;
}

PlacementFrontier::PlacementFrontier(const Circuit& _circ) : circ(_circ) {
  VertexVec input_slice;
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  boolean_in_edges = std::make_shared<b_frontier_t>();

  for (const Qubit& qb : circ.all_qubits()) {
    Vertex input = circ.get_in(qb);
    input_slice.push_back(input);
    Edge candidate = circ.get_nth_out_edge(input, 0);
    quantum_in_edges->insert({qb, circ.skip_irrelevant_edges(candidate)});
  }
  for (const Bit& bit : circ.all_bits()) {
    Vertex input = circ.get_in(bit);
    EdgeVec candidates = circ.get_nth_b_out_bundle(input, 0);
    boolean_in_edges->insert({bit, candidates});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, boolean_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}

void PlacementFrontier::next_slicefrontier() {
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  boolean_in_edges = std::make_shared<b_frontier_t>();
  for (const std::pair<UnitID, Edge>& pair : quantum_out_edges->get<TagKey>()) {
    Edge new_e = circ.skip_irrelevant_edges(pair.second);
    quantum_in_edges->insert({pair.first, new_e});
    Vertex targ = circ.target(new_e);
    EdgeVec targ_classical_ins =
        circ.get_in_edges_of_type(targ, EdgeType::Boolean);
    boolean_in_edges->insert(
        {Bit("frontier_bit", pair.first.index()), targ_classical_ins});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, boolean_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}

QubitGraph monomorph_interaction_graph(
    const Circuit& circ, const unsigned max_edges, unsigned depth_limit) {
  std::set<Qubit> qubits_considered = interacting_qbs(circ);

  QubitGraph q_graph(circ.all_qubits());

  PlacementFrontier current_sf(circ);
  unsigned count_edges = 0;
  for (unsigned slice = 0;
       slice < depth_limit && count_edges < max_edges &&
       !current_sf.slice->empty() && qubits_considered.size() > 1;
       slice++) {
    for (const Vertex& vert : *current_sf.slice) {
      EdgeVec q_out_edges = circ.get_out_edges_of_type(vert, EdgeType::Quantum);
      Qubit qb1;
      Qubit qb2;
      for (const std::pair<UnitID, Edge>& pair :
           current_sf.quantum_out_edges->get<TagKey>()) {
        if (pair.second == q_out_edges[0])
          qb1 = Qubit(pair.first);
        else if (pair.second == q_out_edges[1])
          qb2 = Qubit(pair.first);
      }
      if (!q_graph.edge_exists(qb1, qb2) && !q_graph.edge_exists(qb2, qb1)) {
        q_graph.add_connection(qb1, qb2, slice + 1);
        count_edges++;
      }
    }

    current_sf.next_slicefrontier();
  }
  q_graph.remove_stray_nodes();
  return q_graph;
}

QubitGraph generate_interaction_graph(
    const Circuit& circ, unsigned depth_limit) {
  std::set<Qubit> qubits_considered = interacting_qbs(circ);
  QubitGraph q_graph(circ.all_qubits());
  PlacementFrontier current_sf(circ);

  for (unsigned slice = 0; slice < depth_limit && !current_sf.slice->empty() &&
                           qubits_considered.size() > 1;
       slice++) {
    for (const Vertex& vert : *current_sf.slice) {
      EdgeVec q_out_edges = circ.get_out_edges_of_type(vert, EdgeType::Quantum);
      Qubit qb1;
      Qubit qb2;
      for (const std::pair<UnitID, Edge>& pair :
           current_sf.quantum_out_edges->get<TagKey>()) {
        if (pair.second == q_out_edges[0])
          qb1 = Qubit(pair.first);
        else if (pair.second == q_out_edges[1])
          qb2 = Qubit(pair.first);
      }
      const bool qb1_considered =
          (qubits_considered.find(qb1) != qubits_considered.end());
      const bool qb2_considered =
          (qubits_considered.find(qb2) != qubits_considered.end());
      if ((qb2 != qb1) && (qb1_considered || qb2_considered)) {
        if (!qb1_considered) {
          qubits_considered.erase(qb2);
        } else if (!qb2_considered) {
          qubits_considered.erase(qb1);
        } else if (!q_graph.edge_exists(qb1, qb2)) {
          const unsigned out1 = q_graph.get_degree(qb1);
          const unsigned out2 = q_graph.get_degree(qb2);
          q_graph.add_connection(qb1, qb2, slice + 1);
          if (out1 == 1) {
            qubits_considered.erase(qb1);
          }
          if (out2 == 1) {
            qubits_considered.erase(qb2);
          }
        } else {
          qubits_considered.erase(qb1);
          qubits_considered.erase(qb2);
        }
      }
    }
    current_sf.next_slicefrontier();
  }
  q_graph.remove_stray_nodes();

  return q_graph;
}

QubitLineList qubit_lines(const Circuit& circ) {
  const QubitGraph q_graph = generate_interaction_graph(circ);
  std::set<Qubit> all_qb;
  for (const Qubit& qb : circ.all_qubits()) {
    all_qb.insert(qb);
  }
  QubitLineList found_lines;
  unsigned found_line_size = 0;
  QubitGraph::UndirectedConnGraph graph = q_graph.get_undirected_connectivity();
  do {
    auto u_line = graphs::longest_simple_path(graph);
    QubitLine found;
    std::transform(
        u_line.begin(), u_line.end(), std::back_inserter(found),
        [&graph](auto v) -> Qubit { return graph[v]; });
    found_line_size = found.size();
    if (found_line_size > 1) {
      found_lines.push_back(found);
      for (const auto& vertex : u_line) {
        boost::clear_vertex(vertex, graph);
      }
      for (Qubit q : found) {
        all_qb.erase(q);
      }
    }
  } while (found_line_size > 1);

  for (const Qubit& qb : circ.all_qubits()) {
    if (all_qb.find(qb) != all_qb.end()) found_lines.push_back({qb});
  }
  // print_qubitlines(found_lines);
  return found_lines;
}

// remove a given number of nodes from the architecture, return a set of usable
// nodes (the remainder)
node_set_t best_nodes(Architecture& arc, unsigned n_remove) {
  node_set_t all_nodes = arc.nodes();
  node_set_t bad_nodes;
  // if there are nodes 'removed' already count them as bad nodes
  for (Node n : all_nodes) {
    if (arc.get_degree(n) == 0) {
      bad_nodes.insert(n);
      n_remove--;
    }
  }
  node_set_t removed_nodes = arc.remove_worst_nodes(n_remove);
  bad_nodes.insert(removed_nodes.begin(), removed_nodes.end());
  node_set_t good_nodes;  // keep track of nodes of architecture actually used
  std::set_difference(
      all_nodes.begin(), all_nodes.end(), bad_nodes.begin(), bad_nodes.end(),
      std::inserter(good_nodes, good_nodes.begin()));
  return good_nodes;
}

// map qubit lines to node lines, erasing qubit lines as it goes
qubit_mapping_t map_lines(
    QubitLineList& qb_lines, const std::vector<node_vector_t>& node_lines) {
  qubit_mapping_t outmap;

  // go through all node lines
  for (unsigned i = 0; i < node_lines.size(); i++) {
    node_vector_t node_lst = node_lines[i];
    for (const Node& node : node_lst) {
      outmap.insert({*qb_lines[i].begin(), node});
      qb_lines[i].erase(qb_lines[i].begin());
    }
  }
  return outmap;
}

// trivially place qubit lines on to available nodes, return map
qubit_mapping_t place_qubit_lines(
    const QubitLineList& qb_lines, node_set_t available_nodes) {
  qubit_mapping_t outmap;

  node_set_t::iterator node_it = available_nodes.begin();
  for (const QubitLine& line : qb_lines) {
    for (const Qubit& qb : line) {
      if (node_it == available_nodes.end()) {
        throw ArchitectureInvalidity("Not enough nodes to place all qubits.");
      }
      outmap.insert({qb, *node_it});
      ++node_it;
    }
  }
  return outmap;
}

// Determine intial qubit placement by requesting lines of architecture to place
// lines of qubits on
qubit_mapping_t lines_on_arc(
    Architecture arc, QubitLineList qb_lines, unsigned n_qubits) {
  unsigned difference = arc.n_nodes() - n_qubits;
  // sort from longest to shortest
  std::sort(qb_lines.begin(), qb_lines.end(), [](QubitLine x, QubitLine y) {
    return (x.size() > y.size());
  });

  // get rid of one qubit lines
  while (!qb_lines.empty() && qb_lines.back().size() < 2) {
    difference++;
    qb_lines.pop_back();
  }

  // remove poorly connected nodes, up to the number not used by mapping
  node_set_t unused_nodes = best_nodes(arc, difference);

  // find lengths required
  std::vector<unsigned> lengths;
  for (const QubitLine& line : qb_lines) {
    lengths.push_back(line.size());
  }

  // attempt to find lines of required length on architecture
  std::vector<node_vector_t> node_lines = arc.get_lines(lengths);

  // map qubit lines to node lines to some extent
  qubit_mapping_t outmap = map_lines(qb_lines, node_lines);
  for (qubit_mapping_t::iterator it = outmap.begin(); it != outmap.end();
       ++it) {
    unused_nodes.erase(it->second);
  }
  // map remain remaining qubit lines trivially
  qubit_mapping_t remainder_map = place_qubit_lines(qb_lines, unused_nodes);
  outmap.insert(remainder_map.begin(), remainder_map.end());

  return outmap;
}

// Uses line placement method to produce a suitable qubit_mapping_t for initial
// placement Note that arc is passed by value, since this function modifies it.
qubit_mapping_t line_placement(const Circuit& circ, Architecture arc) {
  QubitLineList qb_lines = qubit_lines(circ);
  if (qb_lines.empty()) {
    return qubit_mapping_t();
  } else {
    unsigned n_qb = circ.n_qubits();
    return lines_on_arc(arc, qb_lines, n_qb);
  }
}

// calculate a cost value for map being considered
double Monomorpher::map_cost(const qubit_bimap_t& n_map) {
  double cost = 0.0;
  const int approx_depth = circ.n_gates() / circ.n_qubits() + 1;
  // constants for scaling single qubit error
  constexpr double c1 = 0.5;
  constexpr double d1 = 1 - 1 / c1;
  for (auto [qb, node] : n_map) {
    // add fidelities of edges from node, weighted by if edge is used by
    // interaction graph
    node_set_t neighs = arc.get_neighbour_nodes(node);
    double edge_sum = 1.0;
    for (Node nei : neighs) {
      double fwd_edge_weighting = 1.0, bck_edge_weighting = 1.0;
      // check if neighbour node is mapped
      r_const_iterator_t nei_qb_it = n_map.right.find(nei);
      if (nei_qb_it == n_map.right.end()) continue;
      // if ( nei_qb_it != n_map.end()){
      unsigned edge_val = q_graph.get_connection_weight(qb, nei_qb_it->second);
      // check if either directed interaction exists
      auto place_interactions_boost = [&](unsigned edge_v) {
        return config.depth_limit - edge_v + 1;
      };
      // if edge is used by interaction in mapping, weight edge higher
      if (edge_val) {
        fwd_edge_weighting += place_interactions_boost(edge_val);
      } else {
        edge_val = q_graph.get_connection_weight(nei_qb_it->second, qb);
        if (edge_val) {
          bck_edge_weighting += place_interactions_boost(edge_val);
        }
      }
      // }
      gate_error_t fwd_error = characterisation.get_error({node, nei});
      gate_error_t bck_error = characterisation.get_error({nei, node});
      edge_sum += fwd_edge_weighting * (1.0 - fwd_error);
      edge_sum += bck_edge_weighting * (1.0 - bck_error);
    }

    // bigger edge sum -> smaller cost
    cost += 1.0 / (edge_sum);

    // add error rate of node
    gate_error_t single_error = characterisation.get_error(node);
    cost += d1 + 1.0 / ((1.0 - single_error) + c1);
    readout_error_t readout_error = characterisation.get_readout_error(node);
    cost += (d1 + 1.0 / ((1.0 - readout_error) + c1)) / (approx_depth * 20);
    // TODO add readout weighting to PlacementConfig?
  }

  return cost;
}

std::vector<MapCost> Monomorpher::place(unsigned max_return) {
  if (max_return < 1)
    throw PlacementError("Max return maps for place must be at least 1.");
  const unsigned interacting_nodes = q_graph.n_connected();
  if ((circ.n_gates() / circ.n_qubits()) >= config.arc_contraction_ratio &&
      circ.n_qubits() > 3) {
    best_nodes(arc, arc.n_nodes() - interacting_nodes);
  }

  std::vector<qubit_bimap_t> potential_maps = monomorphism_edge_break(
      arc, q_graph, config.vf2_max_matches, config.timeout);
  std::vector<MapCost> map_costs;
  for (unsigned i = 0; i < potential_maps.size(); i++) {
    qubit_bimap_t& chosen = potential_maps[i];
    double cost = map_cost(chosen);
    qubit_mapping_t converted_map;
    for (auto [qb, node] : chosen.left) {
      converted_map[qb] = node;
    }
    map_costs.push_back({converted_map, cost});
  }

  max_return = std::min(static_cast<unsigned>(map_costs.size()), max_return);
  std::nth_element(
      map_costs.begin(), map_costs.begin() + max_return - 1, map_costs.end());
  const MapCost nth_map = *(map_costs.begin() + max_return - 1);
  for (std::vector<MapCost>::iterator it = map_costs.begin();
       it != map_costs.end();) {
    if (*it > nth_map) {
      it = map_costs.erase(it);
    } else {
      ++it;
    }
  }

  if (max_return < map_costs.size()) map_costs.resize(max_return);
  return map_costs;
}

}  // namespace tket

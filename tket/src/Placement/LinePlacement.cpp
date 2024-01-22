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

#include "tket/Placement/Placement.hpp"
#include "tket/Utils/HelperFunctions.hpp"

namespace tket {

LinePlacement::LinePlacement(
    const Architecture& _architecture, unsigned _maximum_pattern_gates,
    unsigned _maximum_pattern_depth)
    : GraphPlacement(
          _architecture, 0, 0, _maximum_pattern_gates, _maximum_pattern_depth) {
  architecture_ = _architecture;
}

std::vector<qubit_vector_t> LinePlacement::interactions_to_lines(
    const Circuit& circ_) const {
  qubit_vector_t all_qubits_v = circ_.all_qubits();
  std::vector<WeightedEdge> pattern_weighting =
      this->default_pattern_weighting(circ_);
  if (pattern_weighting.empty()) {
    return {{}};
  }
  const QubitGraph q_graph =
      this->construct_pattern_graph(pattern_weighting, all_qubits_v.size());
  QubitGraph::UndirectedConnGraph uc_graph =
      q_graph.get_undirected_connectivity();

  std::set<Qubit> all_qubits_s(all_qubits_v.begin(), all_qubits_v.end());

  std::vector<qubit_vector_t> found_lines;
  unsigned found_line_size = 0;
  do {
    auto u_line = graphs::longest_simple_path(uc_graph);
    qubit_vector_t found;
    std::transform(
        u_line.begin(), u_line.end(), std::back_inserter(found),
        [&uc_graph](auto v) -> Qubit { return uc_graph[v]; });
    found_line_size = found.size();
    if (found_line_size > 1) {
      found_lines.push_back(found);
      for (const auto& vertex : u_line) {
        boost::clear_vertex(vertex, uc_graph);
      }
      for (const Qubit& q : found) {
        all_qubits_s.erase(q);
      }
    }
  } while (found_line_size > 1);

  for (const Qubit& qb : all_qubits_v) {
    if (all_qubits_s.find(qb) != all_qubits_s.end())
      found_lines.push_back({qb});
  }
  return found_lines;
}

std::map<Qubit, Node> LinePlacement::assign_lines_to_target_graph(
    std::vector<qubit_vector_t>& line_pattern, unsigned n_qubits) const {
  unsigned n_unused_nodes = this->architecture_.n_nodes() - n_qubits;
  // Sort lines from longest to shortest
  std::sort(
      line_pattern.begin(), line_pattern.end(),
      [](qubit_vector_t x, qubit_vector_t y) {
        std::size_t xsz = x.size(), ysz = y.size();
        if (xsz > ysz) {
          return true;
        } else if (xsz < ysz) {
          return false;
        } else {
          return std::lexicographical_compare(
              x.begin(), x.end(), y.begin(), y.end());
        }
      });
  // Remove single qubit lines
  while (!line_pattern.empty() && line_pattern.back().size() < 2) {
    n_unused_nodes++;
    line_pattern.pop_back();
  }
  Architecture copy = this->architecture_;
  node_set_t all_architecture_nodes = this->architecture_.nodes();
  node_set_t bad_nodes;
  for (const Node& node : all_architecture_nodes) {
    if (this->architecture_.get_degree(node) == 0) {
      bad_nodes.insert(node);
      n_unused_nodes--;
    }
  }
  node_set_t removed_nodes = copy.remove_worst_nodes(n_unused_nodes);
  bad_nodes.insert(removed_nodes.begin(), removed_nodes.end());
  node_set_t remaining_nodes;
  std::set_difference(
      all_architecture_nodes.begin(), all_architecture_nodes.end(),
      bad_nodes.begin(), bad_nodes.end(),
      std::inserter(remaining_nodes, remaining_nodes.begin()));

  // store lengths of lines to find
  std::vector<unsigned> lengths;
  for (const qubit_vector_t& line : line_pattern) {
    lengths.push_back(line.size());
  }
  // get lines of Node to map Qubit to

  std::vector<node_vector_t> architecture_lines = copy.get_lines(lengths);
  // and then produce a mapping between qubit and node lines
  std::map<Qubit, Node> out_map;
  for (unsigned i = 0; i < architecture_lines.size(); i++) {
    qubit_vector_t qubit_line = line_pattern[i];
    for (const Node& node : architecture_lines[i]) {
      auto it = qubit_line.begin();
      out_map.insert({*it, node});
      qubit_line.erase(it);
    }
  }
  return out_map;
}

std::vector<std::map<Qubit, Node>> LinePlacement::get_all_placement_maps(
    const Circuit& circ_, unsigned /*matches*/) const {
  std::vector<qubit_vector_t> qb_lines = this->interactions_to_lines(circ_);
  if (!qb_lines.empty()) {
    return {this->assign_lines_to_target_graph(qb_lines, circ_.n_qubits())};
  }
  return {{}};
}

}  // namespace tket

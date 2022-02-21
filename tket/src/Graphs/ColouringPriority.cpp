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

#include "ColouringPriority.hpp"

#include <sstream>
#include <stdexcept>

#include "AdjacencyData.hpp"
#include "Utils/Assert.hpp"

using std::map;
using std::set;
using std::size_t;
using std::string;
using std::vector;

namespace tket {
namespace graphs {

static void fill_initial_node_sequence(
    ColouringPriority::Nodes& nodes, const AdjacencyData& adjacency_data,
    const set<size_t>& vertices_in_component,
    const set<size_t>& initial_clique) {
  nodes.reserve(vertices_in_component.size());
  nodes.clear();

  try {
    for (size_t clique_vertex : initial_clique) {
      // GCOVR_EXCL_START
      if (vertices_in_component.count(clique_vertex) == 0) {
        std::stringstream ss;
        ss << "initial clique vertex " << clique_vertex
           << " is not in this component";
        throw std::runtime_error(ss.str());
      }
      // GCOVR_EXCL_STOP
      nodes.emplace_back();
      nodes.back().vertex = clique_vertex;
    }
    // Now, we do a breadth-first traversal of the remaining vertices,
    // adding all vertices only one step away from the current set.
    size_t current_nodes_begin = 0;

    set<size_t> vertices_seen = initial_clique;
    set<size_t> vertices_to_add;

    for (size_t counter = 0; counter < 2 * vertices_in_component.size();
         ++counter) {
      const size_t current_nodes_end = nodes.size();

      for (size_t i = current_nodes_begin; i < current_nodes_end; ++i) {
        const auto& neighbours = adjacency_data.get_neighbours(nodes[i].vertex);
        for (size_t neighbour : neighbours) {
          if (vertices_seen.count(neighbour) == 0) {
            vertices_to_add.insert(neighbour);
          }
        }
      }
      if (vertices_to_add.empty()) {
        break;
      }
      for (size_t neighbour : vertices_to_add) {
        vertices_seen.insert(neighbour);
        nodes.emplace_back();
        nodes.back().vertex = neighbour;
      }
      vertices_to_add.clear();
      current_nodes_begin = current_nodes_end;
    }
    // GCOVR_EXCL_START
    if (nodes.size() != vertices_in_component.size()) {
      throw std::runtime_error(
          "Final size check: number of filled "
          "nodes does not match number of vertices in this component");
    }
    // GCOVR_EXCL_STOP
  } catch (const std::exception& e) {
    // GCOVR_EXCL_START
    TKET_ASSERT(
        AssertMessage()
        << "ColouringPriority: fill_initial_node_sequence: initial"
        << " clique size " << initial_clique.size() << ", "
        << vertices_in_component.size() << " vertices in"
        << " this component (full graph has "
        << adjacency_data.get_number_of_vertices() << " vertices)."
        << " So far, filled " << nodes.size() << " nodes."
        << " Error: " << e.what());
    // GCOVR_EXCL_STOP
  }
}

// Quadratic, but we're not afraid; the main brute force colouring is
// exponential! Assumes that "fill_initial_node_sequence" has just been called.
// Fills in "earlier_neighbour_node_indices".
static void fill_node_dependencies(
    ColouringPriority::Nodes& nodes, const AdjacencyData& adjacency_data) {
  for (size_t node_index = 1; node_index < nodes.size(); ++node_index) {
    auto& this_node = nodes[node_index];

    for (size_t other_index = 0; other_index < node_index; ++other_index) {
      if (adjacency_data.edge_exists(
              this_node.vertex, nodes[other_index].vertex)) {
        this_node.earlier_neighbour_node_indices.emplace_back(other_index);
      }
    }
  }
}

const ColouringPriority::Nodes& ColouringPriority::get_nodes() const {
  return m_nodes;
}

// GCOVR_EXCL_START
// currently used only within a tket assert macro
string ColouringPriority::print_raw_data(bool relabel_to_simplify) const {
  map<size_t, size_t> old_vertex_to_new_vertex;
  if (relabel_to_simplify) {
    for (size_t i = 0; i < m_nodes.size(); ++i) {
      old_vertex_to_new_vertex[m_nodes[i].vertex] = i;
    }
  } else {
    for (const auto& node : m_nodes) {
      const auto v = node.vertex;
      old_vertex_to_new_vertex[v] = v;
    }
  }
  map<size_t, set<size_t>> data;
  for (const auto& node : m_nodes) {
    const auto this_v = old_vertex_to_new_vertex.at(node.vertex);
    const auto& earlier_v_indices = node.earlier_neighbour_node_indices;
    for (size_t i : earlier_v_indices) {
      const auto other_v = old_vertex_to_new_vertex.at(m_nodes[i].vertex);
      data[this_v].insert(other_v);
      data[other_v].insert(this_v);
    }
  }
  // Compress: remove (i,j) edges if i>j.
  {
    vector<size_t> v_to_erase;
    for (auto& entry : data) {
      v_to_erase.clear();
      for (auto v : entry.second) {
        if (v < entry.first) {
          v_to_erase.push_back(v);
        }
      }
      for (auto v : v_to_erase) {
        entry.second.erase(v);
      }
    }
  }
  std::stringstream ss;

  ss << "\nNeighbours:\nconst std::map<std::size_t, std::vector<std::size_t>> "
        "data { ";
  for (const auto& entry : data) {
    if (!entry.second.empty()) {
      ss << "\n    { " << entry.first << ", { ";
      for (auto v : entry.second) {
        ss << v << ", ";
      }
      ss << "} },";
    }
  }
  ss << "\n};\n\n";
  return ss.str();
}
// GCOVR_EXCL_STOP

ColouringPriority::ColouringPriority(
    const AdjacencyData& adjacency_data,
    const set<size_t>& vertices_in_component, const set<size_t>& initial_clique)
    : m_initial_clique(initial_clique) {
  fill_initial_node_sequence(
      m_nodes, adjacency_data, vertices_in_component, initial_clique);

  fill_node_dependencies(m_nodes, adjacency_data);
}

const set<size_t>& ColouringPriority::get_initial_clique() const {
  return m_initial_clique;
}

}  // namespace graphs
}  // namespace tket

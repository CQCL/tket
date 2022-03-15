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

#include "Path.hpp"

#include <stdexcept>

namespace tket {

template <typename GraphP, typename GraphT>
template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
bool vf2_match_add_callback<GraphP, GraphT>::operator()(
    const CorrespondenceMap1To2 &f, const CorrespondenceMap2To1 &) {
  qubit_bimap_t new_node_map;
  BGL_FORALL_VERTICES_T(v, pattern_graph_, GraphP) {
    Qubit qb = pattern_graph_[v];
    Node node = target_graph_[get(f, v)];
    new_node_map.insert({qb, node});
  }
  n_maps_.push_back(new_node_map);
  return (n_maps_.size() < max);
}

namespace aas {

// The idiomatic way to initialise a PathHandler, and assumes the architecture
// is symmetric. The way using a MatrixXb is for internal use. We initialise
// without using the distance matrix from Architecture, as we generate distances
// using Floyd-Warshall anyway.
PathHandler::PathHandler(const Architecture &arch)
    : PathHandler(arch.get_connectivity()) {}

// breaks for devices with n_qubits >= UINT_MAX/2
// For internal use.
PathHandler::PathHandler(const MatrixXb &connectivity) {
  unsigned n = connectivity.rows();
  size = n;

  unsigned approx_infinity = UINT_MAX >> 1;
  if (n >= approx_infinity) throw std::out_of_range("Qubit number too large");
  distance_matrix_ = MatrixXu::Constant(n, n, approx_infinity);
  path_matrix_ = MatrixXu::Constant(n, n, n);  // set all unreachable nodes to n
  connectivity_matrix_ = connectivity;

  // Floyd-Warshall with path reconstruction, see:
  // https://en.wikipedia.org/wiki/Floyd–Warshall_algorithm#Pseudocode_[11]
  for (unsigned i = 0; i != n; ++i) {
    distance_matrix_(i, i) = 0;
    path_matrix_(i, i) = i;
    for (unsigned j = 0; j != n; ++j) {
      if (i == j) continue;
      if (connectivity_matrix_(i, j)) {
        distance_matrix_(i, j) = 1;
        path_matrix_(i, j) = j;
      }
    }
  }
  for (unsigned k = 0; k != n; ++k) {
    for (unsigned i = 0; i != n; ++i) {
      for (unsigned j = 0; j != n; ++j) {
        unsigned threshold = distance_matrix_(i, k) + distance_matrix_(k, j);
        if (distance_matrix_(i, j) > threshold) {
          distance_matrix_(i, j) = threshold;
          path_matrix_(i, j) = path_matrix_(i, k);
        }
      }
    }
  }
}

PathHandler PathHandler::construct_acyclic_handler() const {
  PathHandler acyclic_handler;
  unsigned n = distance_matrix_.rows();
  MatrixXb acyclic_connectivity(n, n);
  std::vector<unsigned> num_neighbours(n, 0);

  for (unsigned i = 0; i != n; ++i) {
    for (unsigned j = 0; j != n; ++j) {
      if (connectivity_matrix_(i, j)) ++num_neighbours[i];
    }
  }

  for (unsigned i = 0; i != n; ++i) {
    for (unsigned j = 0; j != n; ++j) {
      acyclic_connectivity(i, j) = 0;
    }
  }

  // find the centre
  unsigned threshold_max_dist = n, max_dist_row = 0, centre_node = 0;

  for (unsigned i = 0; i != n; ++i) {
    for (unsigned j = 0; j != n; ++j) {
      if (distance_matrix_(i, j) > max_dist_row)
        max_dist_row = distance_matrix_(i, j);
    }
    if (max_dist_row < threshold_max_dist) {
      centre_node = i;
      threshold_max_dist = max_dist_row;
    }
    max_dist_row = 0;
  }

  // build the acyclic graph outwards from the centre
  std::list<unsigned> current_layer_vertices{centre_node};
  std::list<unsigned> next_layer_vertices;
  std::vector<std::pair<unsigned, unsigned>> parents_neighbours(
      n);  // pair(num_neighbours, parent vertex)
  std::vector<bool> vertices_in_tree(
      n, 0);  // track which vertices are in the acyclic graph
  vertices_in_tree[centre_node] = 1;
  std::pair<unsigned, unsigned> empty_pair{};
  while (!current_layer_vertices.empty()) {
    for (unsigned vert : current_layer_vertices) {
      for (unsigned j = 0; j != n; ++j) {
        if ((distance_matrix_(vert, j) == 1) & (!vertices_in_tree[j])) {
          // if first encounter of vertex, add to next layer
          if (parents_neighbours[j] == empty_pair) {
            next_layer_vertices.push_back(j);
            parents_neighbours[j] = {num_neighbours[vert], vert};

          }
          // choose parent with most neighbours
          else if (parents_neighbours[j].first < num_neighbours[vert]) {
            parents_neighbours[j] = {num_neighbours[vert], vert};
          }
        }
      }
    }
    current_layer_vertices.clear();
    // add in edges
    for (unsigned vert : next_layer_vertices) {
      unsigned parent_vertex = parents_neighbours[vert].second;
      acyclic_connectivity(vert, parent_vertex) = 1;
      acyclic_connectivity(parent_vertex, vert) = 1;

      current_layer_vertices.push_back(vert);
      vertices_in_tree[vert] = 1;
      parents_neighbours[vert] = {};
    }
    next_layer_vertices.clear();
  }

  return PathHandler(acyclic_connectivity);
}

std::list<unsigned> PathHandler::find_path(unsigned i, unsigned j) {
  std::list<unsigned> path{i};
  while (i != j) {
    i = path_matrix_(i, j);
    path.push_back(i);
  }
  return path;
}

std::vector<Node> find_hampath(const Architecture &arch, long timeout) {
  using ArchitectureConn = Architecture::UndirectedConnGraph;
  ArchitectureConn undirected_target = arch.get_undirected_connectivity();
  unsigned n_nodes = arch.n_nodes();

  std::vector<std::pair<Node, Node>> line_nodes(n_nodes - 1);
  for (unsigned n = 0; n != n_nodes - 1; ++n) {
    line_nodes[n] = {Node(n), Node(n + 1)};
  }
  Architecture line_arch(line_nodes);
  ArchitectureConn undirected_pattern = line_arch.get_undirected_connectivity();
  std::vector<qubit_bimap_t> all_maps;
  vf2_match_add_callback<ArchitectureConn, ArchitectureConn> callback(
      all_maps, undirected_pattern, undirected_target, 1);
  bool found_monomorphism = boost::vf2_subgraph_mono(
      undirected_pattern, undirected_target, callback, timeout);

  /* Architecture has no hampath, sad. */
  if (!found_monomorphism) return {};

  /* Left: line, Right: input architecture. */
  const qubit_bimap_t &qmap = all_maps[0];
  std::vector<Node> hampath;
  for (l_const_iterator_t it = qmap.left.begin(); it != qmap.left.end(); ++it) {
    hampath.push_back(it->second);
  }
  return hampath;
}

IterationOrder::IterationOrder(const Architecture &arch) {
  std::set<Node> visited_nodes;

  Node no = *arch.nodes().begin();

  iterationorder.push_back(no);
  visited_nodes.insert(no);

  unsigned whilecount = 0;
  while ((visited_nodes.size() < arch.n_nodes()) &&
         (whilecount < arch.n_nodes())) {
    for (auto edge : arch.get_all_edges_vec()) {
      if (visited_nodes.count(edge.first) &&
          (!visited_nodes.count(edge.second))) {
        iterationorder.push_back(edge.second);
        visited_nodes.insert(edge.second);
        edgelist.push_back(edge);
      } else if (
          (!visited_nodes.count(edge.first)) &&
          (visited_nodes.count(edge.second))) {
        iterationorder.push_back(edge.first);
        visited_nodes.insert(edge.first);
        edgelist.push_back(edge);
      }
    }
    ++whilecount;
  }

  if (visited_nodes.size() != arch.n_nodes()) {
    throw std::logic_error("Unconnected architecture");
  }

  std::reverse(iterationorder.begin(), iterationorder.end());
}

MatrixXb PathHandler::get_connectivity_matrix() const {
  return connectivity_matrix_;
}

MatrixXu PathHandler::get_distance_matrix() const { return distance_matrix_; }

MatrixXu PathHandler::get_path_matrix() const { return path_matrix_; }

unsigned PathHandler::get_size() const { return size; }

}  // namespace aas
}  // namespace tket

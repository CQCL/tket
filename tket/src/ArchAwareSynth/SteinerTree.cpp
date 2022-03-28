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

#include "SteinerTree.hpp"

namespace tket {
namespace aas {

// Steiner Tree generation follows the method of Takahashi and
// Matsuyama et al. see
// https://en.wikipedia.org/wiki/Steiner_tree_problem#Approximating_the_Steiner_tree
SteinerTree::SteinerTree(
    const PathHandler& pathhandler, std::list<unsigned>& nodes_to_add,
    unsigned root_node) {
  root = root_node;
  last_operation_cost = 0;
  init_tree(pathhandler, nodes_to_add);
  while (!nodes_to_add.empty()) {
    add_closest_node_to_tree(pathhandler, nodes_to_add);
  }
  tree_cost = calculate_cost();
}

unsigned SteinerTree::calculate_cost() const {
  unsigned cost = 0;
  for (SteinerNodeType nt : node_types) {
    switch (nt) {
      case SteinerNodeType::ZeroInTree: {
        cost += 2;
        break;
      }
      case SteinerNodeType::OneInTree: {
        ++cost;
        break;
      }
      case SteinerNodeType::Leaf: {
        ++cost;
        break;
      }
      default: {
        // no action occurs
        break;
      }
    }
  }
  if (cost > 0) --cost;
  return cost;
}

unsigned SteinerTree::get_max_element() const {
  unsigned maxelement = 0;
  for (unsigned i = 0; i < node_types.size(); ++i) {
    if (node_types[i] != SteinerNodeType::OutOfTree) maxelement = i;
  }
  return maxelement;
}

std::vector<unsigned> SteinerTree::nodes() const {
  std::vector<unsigned> tree_nodes;
  for (unsigned i = 0; i < node_types.size(); ++i) {
    if (node_types[i] != SteinerNodeType::OutOfTree) {
      tree_nodes.push_back(i);
    }
  }
  return tree_nodes;
}

void SteinerTree::init_tree(
    const PathHandler& pathhandler, std::list<unsigned>& nodes_to_add) {
  if (nodes_to_add.empty())
    throw std::logic_error("Cannot initialise empty Steiner Tree.");

  unsigned n = pathhandler.get_connectivity_matrix().rows();
  node_types = std::vector<SteinerNodeType>(n, SteinerNodeType::OutOfTree);
  num_neighbours = std::vector<unsigned>(n, 0);
  if (nodes_to_add.size() == 1) {
    node_types[nodes_to_add.front()] =
        SteinerNodeType::Leaf;  // A single root is considered a leaf
    tree_nodes = nodes_to_add;
    nodes_to_add.clear();
  } else if (nodes_to_add.size() > 1) {
    // determine the shortest path
    unsigned node1 = nodes_to_add.front(), node2 = nodes_to_add.back();
    unsigned min_distance = pathhandler.get_distance_matrix()(node1, node2);
    for (unsigned no1 : nodes_to_add) {
      for (unsigned no2 : nodes_to_add) {
        if (no1 != no2) {
          unsigned distance = pathhandler.get_distance_matrix()(no1, no2);
          if (distance < min_distance) {
            node1 = no1;
            node2 = no2;
            min_distance = distance;
          }
        }
      }
    }
    // neighbours must both be leaves
    if (pathhandler.get_distance_matrix()(node1, node2) == 1) {
      node_types[node1] = SteinerNodeType::Leaf;
      node_types[node2] = SteinerNodeType::Leaf;
      num_neighbours[node1] = 1;
      num_neighbours[node2] = 1;
      tree_nodes.push_back(node1);
      tree_nodes.push_back(node2);

    }
    // otherwise, 0 nodes exist between.
    else {
      node_types[node1] = SteinerNodeType::Leaf;
      num_neighbours[node1] = 1;
      tree_nodes.push_back(node1);
      add_path_to_tree(pathhandler, node1, node2);
    }
    nodes_to_add.remove(node1);
    nodes_to_add.remove(node2);
  }
}

void SteinerTree::add_path_to_tree(
    const PathHandler& pathhandler, unsigned node_in_tree,
    unsigned node_to_add) {
  // last node is a leaf, with one neighbor
  node_types[node_to_add] = SteinerNodeType::Leaf;
  num_neighbours[node_to_add] = 1;
  tree_nodes.push_back(node_to_add);

  if ((node_in_tree == pathhandler.get_size()) ||
      (node_to_add == pathhandler.get_size())) {
    throw std::logic_error("searching for a node which is not in the tree");
  }

  if (pathhandler.get_distance_matrix()(node_in_tree, node_to_add) <
      pathhandler.get_distance_matrix()(node_to_add, node_in_tree)) {
    if (pathhandler.get_path_matrix()(node_to_add, node_in_tree) ==
        pathhandler.get_size()) {
      node_in_tree = pathhandler.get_path_matrix()(node_in_tree, node_to_add);
    } else {
      node_in_tree = pathhandler.get_path_matrix()(node_to_add, node_in_tree);
    }

    if ((node_in_tree == pathhandler.get_size()) ||
        (node_to_add == pathhandler.get_size())) {
      throw std::logic_error("searching for a node which is not in the tree");
    }

    while (node_in_tree != node_to_add) {
      // node in a path = zero node, it has two neighbbors
      node_types[node_in_tree] = SteinerNodeType::ZeroInTree;
      tree_nodes.push_back(node_in_tree);
      num_neighbours[node_in_tree] = 2;

      if (pathhandler.get_path_matrix()(node_to_add, node_in_tree) ==
          pathhandler.get_size()) {
        node_in_tree = pathhandler.get_path_matrix()(node_in_tree, node_to_add);
      } else {
        node_in_tree = pathhandler.get_path_matrix()(node_to_add, node_in_tree);
      }

      if ((node_in_tree == pathhandler.get_size()) ||
          (node_to_add == pathhandler.get_size())) {
        throw std::logic_error("searching for a node which is not in the tree");
      }
    }
  } else {
    if (pathhandler.get_path_matrix()(node_to_add, node_in_tree) ==
        pathhandler.get_size()) {
      node_to_add = pathhandler.get_path_matrix()(node_in_tree, node_to_add);
    } else {
      node_to_add = pathhandler.get_path_matrix()(node_to_add, node_in_tree);
    }

    if ((node_in_tree == pathhandler.get_size()) ||
        (node_to_add == pathhandler.get_size())) {
      throw std::logic_error("searching for a node which is not in the tree");
    }

    while (node_in_tree != node_to_add) {
      node_types[node_to_add] =
          SteinerNodeType::ZeroInTree;  // node in a path = zero node

      tree_nodes.push_back(node_to_add);

      num_neighbours[node_to_add] = 2;  // it has two neighbbors

      if (pathhandler.get_path_matrix()(node_to_add, node_in_tree) ==
          pathhandler.get_size()) {
        node_to_add = pathhandler.get_path_matrix()(node_in_tree, node_to_add);
      } else {
        node_to_add = pathhandler.get_path_matrix()(node_to_add, node_in_tree);
      }

      if ((node_in_tree == pathhandler.get_size()) ||
          (node_to_add == pathhandler.get_size())) {
        throw std::logic_error("searching for a node which is not in the tree");
      }
    }
  }
}

void SteinerTree::add_closest_node_to_tree(
    const PathHandler& pathhandler, std::list<unsigned>& nodes_to_add) {
  unsigned closest_node = 0;
  unsigned distant_node = UINT_MAX;  // initialise to "infinity"
  unsigned closest_tree_node(
      tree_nodes.front());  // closest steiner point to the nodes we want to add

  std::list<unsigned>::iterator itr1;
  std::list<unsigned>::iterator itr2;

  // take the two closest points (one of the tree, the other among the
  // nodes_to_add)
  for (unsigned node_to_add : nodes_to_add) {
    for (unsigned tree_node : tree_nodes) {
      if (pathhandler.get_distance_matrix()(tree_node, node_to_add) <
          distant_node)  // graph is oriented, we want distance
                         // tree->nodes_to_add
      {
        closest_node = node_to_add;
        closest_tree_node = tree_node;
        distant_node =
            pathhandler.get_distance_matrix()(tree_node, node_to_add);
      }
    }
  }
  nodes_to_add.remove(closest_node);

  // if the node of the tree is a leaf, then its a one node after adding the
  // new node
  if (node_types[closest_tree_node] == SteinerNodeType::Leaf) {
    node_types[closest_tree_node] = SteinerNodeType::OneInTree;
  }
  ++num_neighbours[closest_tree_node];
  add_path_to_tree(pathhandler, closest_tree_node, closest_node);
}

// assumes i & j are neighbouring vertices
int SteinerTree::cost_of_operation(unsigned i, unsigned j) const {
  SteinerNodeType i_type = node_types[i];
  SteinerNodeType j_type = node_types[j];

  // tedious case analysis
  switch (i_type) {
    case SteinerNodeType::ZeroInTree: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree:
        case SteinerNodeType::OneInTree:
        case SteinerNodeType::Leaf:
        case SteinerNodeType::OutOfTree: {
          return 0;
        }
        default: {
          throw InvalidCostCalculation(
              "[AAS]: Invalid cost calculation, wrong SteinerNodeType");
        }
      }
    }
    case SteinerNodeType::OneInTree: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree:
        case SteinerNodeType::Leaf: {
          return -1;
        }
        case SteinerNodeType::OneInTree:
        case SteinerNodeType::OutOfTree: {
          return 1;
        }
        default: {
          throw InvalidCostCalculation(
              "[AAS]: Invalid cost calculation, wrong SteinerNodeType");
        }
      }
    }
    case SteinerNodeType::Leaf: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree:
        case SteinerNodeType::Leaf: {
          return -1;
        }
        case SteinerNodeType::OneInTree:
        case SteinerNodeType::OutOfTree: {
          return 1;
        }
        default: {
          throw InvalidCostCalculation(
              "[AAS]: Invalid cost calculation, wrong SteinerNodeType");
        }
      }
    }
    case SteinerNodeType::OutOfTree: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree:
        case SteinerNodeType::OneInTree:
        case SteinerNodeType::Leaf:
        case SteinerNodeType::OutOfTree: {
          return 0;
        }
        default: {
          throw InvalidCostCalculation(
              "[AAS]: Invalid cost calculation, wrong SteinerNodeType");
        }
      }
    }
    default: {
      throw InvalidCostCalculation(
          "[AAS]: Invalid cost calculation, wrong SteinerNodeType");
    }
  }
}

OperationList SteinerTree::operations_available(
    const PathHandler& pathhandler) const {
  OperationList operations;
  for (unsigned i = 0; i != node_types.size(); ++i) {
    for (unsigned j = 0; j != node_types.size(); ++j) {
      if (i == j) continue;
      if (!pathhandler.get_connectivity_matrix()(i, j)) continue;

      if (node_types[i] == SteinerNodeType::OneInTree ||
          node_types[i] == SteinerNodeType::Leaf) {
        if (node_types[j] == SteinerNodeType::ZeroInTree ||
            node_types[j] == SteinerNodeType::Leaf) {
          operations.push_back({i, j});
        }
      }
    }
  }
  return operations;
}

void SteinerTree::add_row(unsigned i, unsigned j) {
  /* CNOT with control j and target i */
  SteinerNodeType i_type = node_types[i];
  SteinerNodeType j_type = node_types[j];

  last_operation_cost = cost_of_operation(i, j);
  tree_cost += last_operation_cost;

  switch (i_type) {
    case SteinerNodeType::ZeroInTree: {
      // no action occurs
      break;
    }
    case SteinerNodeType::OneInTree: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree: {
          node_types[j] = SteinerNodeType::OneInTree;
          break;
        }
        case SteinerNodeType::OneInTree: {
          node_types[j] = SteinerNodeType::ZeroInTree;
          break;
        }
        case SteinerNodeType::Leaf: {
          TKET_ASSERT(num_neighbours[i] != 0);
          TKET_ASSERT(num_neighbours[j] != 0);
          node_types[j] = SteinerNodeType::OutOfTree;
          num_neighbours[i] -= 1;
          num_neighbours[j] -= 1;
          if (num_neighbours[i] == 1) {
            node_types[i] = SteinerNodeType::Leaf;
          }
          break;
        }
        case SteinerNodeType::OutOfTree: {
          node_types[j] = SteinerNodeType::Leaf;
          node_types[i] = SteinerNodeType::OneInTree;
          num_neighbours[i] += 1;
          num_neighbours[j] += 1;
          break;
        }
        default: {
          throw InvalidRowOperation(
              "[AAS]: Invalid row operation, invalid combination "
              "SteinerNodeType in add_row");
        }
      }
      break;
    }
    case SteinerNodeType::Leaf: {
      switch (j_type) {
        case SteinerNodeType::ZeroInTree: {
          node_types[j] = SteinerNodeType::OneInTree;
          break;
        }
        case SteinerNodeType::OneInTree: {
          node_types[j] = SteinerNodeType::ZeroInTree;
          break;
        }
        // only happens when there are 2 vertices lefts
        case SteinerNodeType::Leaf: {
          TKET_ASSERT(num_neighbours[i] != 0);
          TKET_ASSERT(num_neighbours[j] != 0);
          node_types[j] = SteinerNodeType::OutOfTree;
          node_types[i] = SteinerNodeType::OutOfTree;
          num_neighbours[i] -= 1;
          num_neighbours[j] -= 1;
          break;
        }
        case SteinerNodeType::OutOfTree: {
          node_types[j] = SteinerNodeType::Leaf;
          node_types[i] = SteinerNodeType::OneInTree;
          num_neighbours[i] += 1;
          num_neighbours[j] += 1;
          break;
        }
        default: {
          // Invalid combination of nodes types in add row operation
          TKET_ASSERT(false);
        }
      }
    }
    case SteinerNodeType::OutOfTree: {
      // no action occurs
      break;
    }
    default: {
      TKET_ASSERT(!"Invalid combination of nodes types in add row operation");
    }
  }
}

bool SteinerTree::fully_reduced() const { return tree_cost == 0; }

static std::pair<unsigned, std::vector<unsigned>> steiner_reduce(
    Circuit& circ, DiagMatrix& CNOT_matrix, const PathHandler& paths,
    unsigned col, unsigned root, std::list<unsigned>& nodes, bool upper,
    CNotSynthType cnottype) {
  std::pair<unsigned, std::vector<unsigned>> result;

  PathHandler directed_paths;

  std::list<unsigned> fresh_node_list(nodes);

  if (upper) {
    MatrixXb directed_connectivity = paths.get_connectivity_matrix();
    for (unsigned i = 0; i < directed_connectivity.rows(); ++i) {
      for (unsigned j = 0; j < directed_connectivity.cols(); ++j) {
        if (i < root) directed_connectivity(i, j) = 0;
        if (j < root) directed_connectivity(i, j) = 0;
      }
    }

    directed_paths = PathHandler(directed_connectivity);
  } else {
    MatrixXb directed_connectivity = paths.get_connectivity_matrix();
    if (cnottype == CNotSynthType::HamPath) {
      for (unsigned i = 0; i < directed_connectivity.rows(); ++i) {
        for (unsigned j = 0; j < directed_connectivity.cols(); ++j) {
          if (j > i) directed_connectivity(i, j) = 0;
          if ((j + 1 != i) && (j != i + 1)) directed_connectivity(i, j) = 0;
        }
      }
    }

    directed_paths = PathHandler(directed_connectivity);
  }
  SteinerTree cnot_tree = SteinerTree(directed_paths, fresh_node_list, root);

  // make list of edges, starting from root to leaves of tree
  std::list<std::pair<unsigned, unsigned>> parent_child_list;
  std::set<unsigned> possible_parents{root};
  unsigned n_edges = cnot_tree.tree_nodes.size();
  if (n_edges > 0) --n_edges;
  std::set<unsigned> visited_parents{root};
  unsigned counterwhile = 0;
  while ((parent_child_list.size() < n_edges) &&
         (counterwhile < n_edges * n_edges)) {
    ++counterwhile;
    std::set<unsigned> new_parents;

    for (unsigned node : cnot_tree.tree_nodes) {
      for (unsigned parent : possible_parents) {
        if (directed_paths.get_connectivity_matrix()(parent, node)) {
          if (visited_parents.find(node) == visited_parents.end()) {
            new_parents.insert(node);
            visited_parents.insert(node);
            parent_child_list.push_back({parent, node});
          }
        }
      }
    }

    possible_parents = new_parents;
  }

  /* Remove zeros */
  if (upper) {
    std::list<std::pair<unsigned, unsigned>> zeros;
    for (auto& [s0, s1] : parent_child_list) {
      if (!CNOT_matrix._matrix(s0, col)) {
        zeros.push_back({s0, s1});
      }
    }
    while (!zeros.empty()) {
      std::pair<unsigned, unsigned> s_pair = zeros.back();
      unsigned s0 = s_pair.first, s1 = s_pair.second;
      zeros.pop_back();
      if (!CNOT_matrix._matrix(s0, col)) {
        CNOT_matrix.row_add(s1, s0);
        circ.add_op<unsigned>(OpType::CX, {s1, s0});
      }
    }
    parent_child_list.reverse();
    for (auto& [s0, s1] : parent_child_list) {
      CNOT_matrix.row_add(s0, s1);
      circ.add_op<unsigned>(OpType::CX, {s0, s1});
    }
  } else {
    for (auto& [s0, s1] : parent_child_list) {
      if (!CNOT_matrix._matrix(s1, col)) {
        CNOT_matrix.row_add(s0, s1);
        circ.add_op<unsigned>(OpType::CX, {s0, s1});
      }
    }
    parent_child_list.reverse();
    for (auto& [s0, s1] : parent_child_list) {
      CNOT_matrix.row_add(s0, s1);
      circ.add_op<unsigned>(OpType::CX, {s0, s1});
    }
  }
  result.first = cnot_tree.get_max_element();
  result.second = cnot_tree.nodes();
  return result;
}

static std::pair<unsigned, std::vector<unsigned>> steiner_reduce_rec(
    Circuit& circ, DiagMatrix& CNOT_matrix, const PathHandler& paths,
    unsigned col, unsigned root, std::list<unsigned>& nodes) {
  std::pair<unsigned, std::vector<unsigned>> result;

  std::list<unsigned> fresh_node_list(nodes);

  SteinerTree cnot_tree = SteinerTree(paths, fresh_node_list, root);

  // make list of edges, starting from root to leaves of tree
  std::list<std::pair<unsigned, unsigned>> parent_child_list;
  std::set<unsigned> possible_parents{root};
  unsigned n_edges = cnot_tree.tree_nodes.size();
  if (n_edges > 0) --n_edges;
  std::set<unsigned> visited_parents{root};
  unsigned counterwhile = 0;
  while ((parent_child_list.size() < n_edges) &&
         (counterwhile < n_edges * n_edges)) {
    ++counterwhile;
    std::set<unsigned> new_parents;

    for (unsigned node : cnot_tree.tree_nodes) {
      for (unsigned parent : possible_parents) {
        if (paths.get_connectivity_matrix()(parent, node)) {
          if (visited_parents.find(node) == visited_parents.end()) {
            new_parents.insert(node);
            visited_parents.insert(node);
            parent_child_list.push_back({parent, node});
          }
        }
      }
    }

    possible_parents = new_parents;
  }

  /* Remove zeros */
  for (auto& [s0, s1] : parent_child_list) {
    if (!CNOT_matrix._matrix(s1, col)) {
      CNOT_matrix.row_add(s0, s1);
      circ.add_op<unsigned>(OpType::CX, {s0, s1});
    }
  }
  parent_child_list.reverse();
  for (auto& [s0, s1] : parent_child_list) {
    CNOT_matrix.row_add(s0, s1);
    circ.add_op<unsigned>(OpType::CX, {s0, s1});
  }
  result.first = cnot_tree.get_max_element();
  result.second = cnot_tree.nodes();
  return result;
}

void aas_cnot_synth_rec(
    DiagMatrix& CNOT_matrix, const PathHandler& paths,
    std::vector<unsigned>& pivot_cols, Circuit& cnot_circuit,
    std::vector<unsigned> usablenodes = {}) {
  // order usable nodes from highest to lowest element, the list should already
  // be close to the opposite of this order anyway.
  std::sort(usablenodes.begin(), usablenodes.end());
  std::reverse(usablenodes.begin(), usablenodes.end());
  for (unsigned current_row : usablenodes) {
    unsigned pivot = pivot_cols[current_row];

    std::list<unsigned> nodes;

    for (unsigned r = 0; r != current_row; ++r) {
      if (CNOT_matrix._matrix(r, pivot)) {
        nodes.push_back(r);
      }
    }
    unsigned max_node_in_tree = 0;
    std::vector<unsigned> new_usable_nodes;
    if (!nodes.empty()) {
      nodes.push_back(current_row);
      std::pair<unsigned, std::vector<unsigned>> res = steiner_reduce_rec(
          cnot_circuit, CNOT_matrix, paths, pivot, current_row, nodes);
      max_node_in_tree = res.first;
      new_usable_nodes = res.second;
    }

    if (max_node_in_tree > current_row) {
      aas_cnot_synth_rec(
          CNOT_matrix, paths, pivot_cols, cnot_circuit, new_usable_nodes);
    }
  }
}

// see https://arxiv.org/abs/2004.06052 and https://github.com/Quantomatic/pyzx
// for more information.
Circuit aas_CNOT_synth(
    DiagMatrix& CNOT_matrix, const PathHandler& paths, CNotSynthType cnottype) {
  unsigned pivot = 0;
  unsigned max_node_in_tree = 0;
  std::vector<unsigned> usable_nodes(paths.get_size());
  std::iota(usable_nodes.begin(), usable_nodes.end(), 0);

  std::vector<unsigned> pivot_cols;
  Circuit cnot_circuit(paths.get_size());
  for (unsigned current_row = 0; current_row != CNOT_matrix.n_rows();
       ++current_row) {
    bool found_pivot = false;
    std::list<unsigned> nodes;
    while ((!found_pivot) && pivot < CNOT_matrix.n_cols()) {
      for (unsigned r = current_row; r != CNOT_matrix.n_rows(); ++r) {
        if (CNOT_matrix._matrix(r, pivot)) {
          nodes.push_back(r);
        }
      }
      if (!nodes.empty()) {
        pivot_cols.push_back(pivot);
        found_pivot = true;
      } else {
        ++pivot;
      }
    }

    // can't try any more pivots
    if (!found_pivot)
      throw std::logic_error("Could not find pivot node in CNOT synthesis.");

    bool row_not_in_notes = std::none_of(
        nodes.cbegin(), nodes.cend(),
        [current_row](unsigned no) { return no == current_row; });

    if (row_not_in_notes) nodes.push_front(current_row);
    std::pair<unsigned, std::vector<unsigned>> res = steiner_reduce(
        cnot_circuit, CNOT_matrix, paths, pivot, current_row, nodes, true,
        cnottype);

    max_node_in_tree = res.first;
    usable_nodes = res.second;
    ++pivot;
  }

  unsigned current_row = pivot_cols.size() - 1;
  while (current_row > 0) {
    std::list<unsigned> nodes;
    if (CNOT_matrix.is_id_until_columns(current_row)) {
      pivot = pivot_cols[current_row];

      for (unsigned r = 0; r != current_row; ++r) {
        if (CNOT_matrix._matrix(r, pivot)) {
          nodes.push_back(r);
        }
      }
      if (!nodes.empty()) {
        nodes.push_back(current_row);

        std::pair<unsigned, std::vector<unsigned>> res = steiner_reduce(
            cnot_circuit, CNOT_matrix, paths, pivot, current_row, nodes, false,
            cnottype);

        max_node_in_tree = res.first;
        usable_nodes = res.second;
      } else {
        max_node_in_tree = 0;
        std::iota(usable_nodes.begin(), usable_nodes.end(), 0);
      }
    }

    if ((max_node_in_tree > current_row) && (cnottype == CNotSynthType::Rec)) {
      aas_cnot_synth_rec(
          CNOT_matrix, paths, pivot_cols, cnot_circuit, usable_nodes);
    }

    --current_row;

    TKET_ASSERT(CNOT_matrix.is_id_until_columns(current_row));
  }

  // when constructing circuit, we want to
  // prepend gates, but it is easier to append and then take the dagger.
  return cnot_circuit;
}

Circuit CNotSwapSynth::get_circuit() { return circ; }

void CNotSwapSynth::add_swap(unsigned first, unsigned second) {
  CNOT_matrix.row_add(first, second);
  CNOT_matrix.row_add(second, first);
  CNOT_matrix.row_add(first, second);
  circ.add_op<unsigned>(OpType::CX, {first, second});
  circ.add_op<unsigned>(OpType::CX, {second, first});
  circ.add_op<unsigned>(OpType::CX, {first, second});
}

Circuit aas_CNOT_synth_SWAP(DiagMatrix& CNOT_matrix, const PathHandler& paths) {
  CNotSwapSynth cnot(paths, CNOT_matrix);
  TKET_ASSERT(cnot.valid_result());
  return cnot.get_circuit();
}

CNotSwapSynth::CNotSwapSynth(
    const PathHandler& pathhandler, const DiagMatrix& CNOT_mat)
    : paths(pathhandler),
      CNOT_matrix(CNOT_mat),
      circ(Circuit(paths.get_size())) {
  for (unsigned current_row = 0; current_row != CNOT_matrix.n_rows();
       ++current_row) {
    if (!CNOT_matrix._matrix(current_row, current_row)) {
      // try to find element equal to 1
      // and set diagonal element to 1
      unsigned one = current_row;
      while (!CNOT_matrix._matrix(one, current_row)) {
        ++one;
      }
      unsigned current_node = swap_to_root(one, current_row);

      // remove the 1 with the use of the root
      CNOT_matrix.row_add(current_node, current_row);
      circ.add_op<unsigned>(OpType::CX, {current_node, current_row});
      cleanup_swaps();
    }

    if (!CNOT_matrix._matrix(current_row, current_row)) {
      throw std::logic_error(
          "The given matrix is not invertible, the input was not created by a "
          "cnot circuit");
    }

    for (unsigned r = current_row + 1; r != CNOT_matrix.n_rows(); ++r) {
      if (CNOT_matrix._matrix(r, current_row)) {
        unsigned current_node = swap_to_root(r, current_row);

        CNOT_matrix.row_add(current_row, current_node);
        circ.add_op<unsigned>(OpType::CX, {current_row, current_node});
        cleanup_swaps();
      }
    }
  }
  // end of upper creation

  for (unsigned current_row = CNOT_matrix.n_rows() - 1; current_row > 0;
       --current_row) {
    for (unsigned r = 0; r < current_row; ++r) {
      if (CNOT_matrix._matrix(r, current_row)) {
        unsigned current_node = swap_to_root(r, current_row);

        CNOT_matrix.row_add(current_row, current_node);
        circ.add_op<unsigned>(OpType::CX, {current_row, current_node});

        cleanup_swaps();
      }
    }
  }
}

void CNotSwapSynth::cleanup_swaps() {
  while (!swaps.empty()) {
    std::pair<unsigned, unsigned> swap = swaps.top();
    swaps.pop();
    add_swap(swap.first, swap.second);
  }
}

unsigned CNotSwapSynth::swap_to_root(
    unsigned start_node, unsigned current_row) {
  unsigned current_node = start_node;

  while (paths.get_path_matrix()(current_node, current_row) != current_row) {
    unsigned new_node = paths.get_path_matrix()(current_node, current_row);
    add_swap(current_node, new_node);

    swaps.push({current_node, new_node});
    current_node = new_node;
  }
  return current_node;
}

bool CNotSwapSynth::valid_result() { return CNOT_matrix.is_id(); }

}  // namespace aas
}  // namespace tket

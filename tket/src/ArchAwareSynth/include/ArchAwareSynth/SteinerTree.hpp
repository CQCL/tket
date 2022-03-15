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

#pragma once
#include "Converters/Gauss.hpp"
#include "Path.hpp"
#include "Utils/Exceptions.hpp"

namespace tket {
namespace aas {

typedef std::pair<unsigned, unsigned> Operation;
typedef std::list<Operation> OperationList;

enum class CNotSynthType { SWAP, HamPath, Rec };

/**
 * Nodes in a Steiner tree may correspond to different values in the biadjacency
 * action matrix or phase polynomial
 */
enum class SteinerNodeType {
  // "in tree" := not a leaf & not disconnected from tree
  ZeroInTree,
  OneInTree,
  // all leaves are 1s by definition
  Leaf,
  OutOfTree
};

class InvalidCostCalculation : public std::logic_error {
 public:
  explicit InvalidCostCalculation(const std::string &message)
      : std::logic_error(message) {}
};

class InvalidRowOperation : public std::logic_error {
 public:
  explicit InvalidRowOperation(const std::string &message)
      : std::logic_error(message) {}
};

/**
 * This clas is creating a steiner tree, which is a mst including all the nodes
 * of a given phase gadget plus all the nodes which are neede to connect all the
 * nodes of the gadget. This class also offers the function to reduce this tree
 * by extracting operations steps by step. This class is designed with
 * architecture aware synthesis in mind; it may not be suitable for other
 * generic purposes.
 */
class SteinerTree {
 public:
  SteinerTree() {}
  /**
   * Construct a Steinertree from given parameters
   * @param pathhandler giving the included edges
   * @param nodes_to_add list by ref, will be changed during the construction of
   * the tree
   * @param root_node initial node, deleted in the last step
   */
  explicit SteinerTree(
      const PathHandler &pathhandler, std::list<unsigned> &nodes_to_add,
      unsigned root_node);  // N.B. non-const list by ref!

  /**
   * calculate cost of performing a CNOT between two neighbouring nodes,
   * dependent on the SteinerNodeTypes.
   * @param i control index of operation
   * @param j target index of operation
   */
  int cost_of_operation(unsigned i, unsigned j) const;

  /**
   * initialises the tree
   * @param pathhandler connectivity for the construction
   * @param nodes_to_add list of nodes which should be added, the list will be
   * changes during the process
   *
   */
  void init_tree(
      const PathHandler &pathhandler, std::list<unsigned> &nodes_to_add);

  /**
   * adds the closesd nodes to a given list to the tree
   * @param pathhandler connectivity for the construction
   * @param nodes_to_add list of nodes which should be added, the list will be
   */
  void add_closest_node_to_tree(
      const PathHandler &pathhandler, std::list<unsigned> &nodes_to_add);

  /**
   * adds a list of nodes to a given list to the tree
   * @param pathhandler connectivity for the construction
   * @param node_in_tree node in tree which is closest to the node to add
   * @param node_to_add node which should be added
   */
  void add_path_to_tree(
      const PathHandler &pathhandler, unsigned node_in_tree,
      unsigned node_to_add);

  /**
   * calculates the available operation in this tree which could be executed
   * @param pathhandler conections used for the calculation
   * @return gives the list of all avilable operations
   */
  OperationList operations_available(const PathHandler &pathhandler) const;

  /**
   * Implements a CNOT from node j to node i, updates costs and tree
   * @param i control index
   * @param j target index
   */
  void add_row(unsigned i, unsigned j);

  /**
   * checks is the tree is fully reduced
   * @return true if fully reduced
   */
  bool fully_reduced() const;

  /**
   * gives the cost of the tree
   * @return cost of the tree
   */
  unsigned calculate_cost() const;

  /**
   * returns the nodes of the tree which has the highest index. only the out of
   * tree nodes are not taken into consideration
   * @return index of maximum node
   */
  unsigned get_max_element() const;

  /**
   * gives all nodes of the tree which are leaf, one in or zero in
   * @return ordered list of all nodes
   */
  std::vector<unsigned> nodes() const;

  unsigned tree_cost;  // the cost to reduce the Steiner tree alone
  int last_operation_cost;
  unsigned root;
  std::vector<SteinerNodeType> node_types;
  std::vector<unsigned> num_neighbours;
  std::list<unsigned> tree_nodes;
};

/**
 * object to store and perform the swap nased cnot synthesis
 */
class CNotSwapSynth {
 public:
  CNotSwapSynth() {}
  /**
   * construct the object for the cnot swap synth and perform the reduction
   * @param pathhandler path for reduction
   * @param CNOT_mat input of executed cnots which should be reduced
   */
  explicit CNotSwapSynth(
      const PathHandler &pathhandler, const DiagMatrix &CNOT_mat);

  /**
   * gives the calculated circuit
   */
  Circuit get_circuit();

  /**
   * checks if the matrix is id after reduction
   * @return true if matrix is id
   */
  bool valid_result();

 private:
  void add_swap(unsigned first, unsigned second);
  void cleanup_swaps();
  unsigned swap_to_root(unsigned start_node, unsigned current_row);

  PathHandler paths;
  DiagMatrix CNOT_matrix;
  Circuit circ;
  std::stack<std::pair<unsigned, unsigned>> swaps;
};

/**
 * This method uses Kissinger & de Meijer's Steiner-Gauss
 * (https://arxiv.org/abs/1904.00633), and reduces the CNOT_matrix, which is
 * non-const. This function offers the recursive algorithm and the Hamilton
 * based algorithm
 * @param CNOT_matrix input of executed cnots which should be reduced
 * @param paths pathhandler used for the reduction
 * @param cnottype type of algorithm which could be CNotSynthType::Rec or
 * CNotSynthType::HamPath
 * @return routed circuit
 */
Circuit aas_CNOT_synth(
    DiagMatrix &CNOT_matrix, const PathHandler &paths,
    CNotSynthType cnottype = CNotSynthType::Rec);

/**
 * this method offers the swap bases cnot synth
 * @param CNOT_matrix input of executed cnots which should be reduced
 * @param paths pathhandler used for the reduction
 * @return routed circuit
 */
Circuit aas_CNOT_synth_SWAP(DiagMatrix &CNOT_matrix, const PathHandler &paths);

}  // namespace aas
}  // namespace tket

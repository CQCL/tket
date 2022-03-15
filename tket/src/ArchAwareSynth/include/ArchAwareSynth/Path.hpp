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
#include "Architecture/Architecture.hpp"
#include "Placement/Placement.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/UnitID.hpp"
namespace tket {
namespace aas {

typedef Eigen::Matrix<unsigned, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    MatrixXu;

/**
 *  Holds distances & paths between nodes -- can optionally remove edges to form
 * a tree.
 */
class PathHandler {
 public:
  PathHandler() {}

  /**
   * Initialises the pathhandler with a given architecture. Architecture
   * initialisation assumes symmetric connectivity.
   * @param arch architecture giving the undirected graph for the initialisation
   * graph
   */
  explicit PathHandler(const Architecture &arch);

  /**
   * Initialises the pathhandler with a given connectivity matrix. This function
   * should only be used in the aas code. This function interprets the matrix as
   * a directed graph.
   * @param connectivity matrix with the different connections of the directed
   * graph
   */
  explicit PathHandler(const MatrixXb &connectivity);

  /**
   * Returns a handler for a spanning tree of the architecture.
   */
  PathHandler construct_acyclic_handler() const;

  /**
   * Find a path between two vertices in the architecture.
   * @param i start node for the path
   * @param j end node for the path
   */
  std::list<unsigned> find_path(unsigned i, unsigned j);

  /**
   * get connectivity_matrix_ of the pathhandler
   * @return connectivity_matrix_
   */
  MatrixXb get_connectivity_matrix() const;

  /**
   * get distance_matrix_ of the pathhandler
   * @return distance_matrix_
   */
  MatrixXu get_distance_matrix() const;

  /**
   * get path_matrix_ of the pathhandler
   * @return path_matrix_
   */
  MatrixXu get_path_matrix() const;

  /**
   * get size of the pathhandler
   * @return size
   */
  unsigned get_size() const;

 private:
  MatrixXb connectivity_matrix_;
  MatrixXu distance_matrix_;
  MatrixXu path_matrix_;
  unsigned size;
};

/**
 * class to calculate and store the iteration order needed for the cnot synth in
 * the architecture aware synth code. The iteration order makes sure, that the
 * not yet iterated nodes are still connected.
 */
class IterationOrder {
 public:
  IterationOrder() {}
  /**
   * construct and calculate the iterationorder.
   * @param arch given architecture
   */
  explicit IterationOrder(const Architecture &arch);

  /**
   * give the in the constructur calculated iteration order
   * @return ordered vector of the architecutre in which they can be iterated.
   */
  std::vector<Node> get_iterationorder() { return iterationorder; }

  /**
   * give the in the constructur calculated edges used for the iteration
   * @return vector of edges pf the architecutre which are used for the
   * iteration.
   */
  std::vector<std::pair<Node, Node>> get_edgelist() { return edgelist; }

 private:
  std::vector<Node> iterationorder;
  std::vector<std::pair<Node, Node>> edgelist;
};

/**
 * Find a Hamiltonian path in the architecture. Returns {} if no
 * Hamiltonian path is found within timeout. Timeout is in ms.
 *
 * @param arch architecture where the path is searched
 * @param timeout give the timeout for the search of the hamiltonpath,
 *                  default value is 10000
 * @return ordered vector of nodes in the path
 */
std::vector<Node> find_hampath(const Architecture &arch, long timeout = 10000);

}  // namespace aas
}  // namespace tket

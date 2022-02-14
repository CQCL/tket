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
#include "Utils/BiMapHeaders.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

typedef std::map<Node, Node> interacting_nodes_t;
typedef std::pair<Node, Node> swap_t;
typedef std::set<swap_t> swap_set_t;
typedef std::vector<size_t> lexicographical_distances_t;

class LexicographicalComparisonError : public std::logic_error {
 public:
  explicit LexicographicalComparisonError(const std::string& message)
      : std::logic_error(message) {}
};

/**
 * A class for running lexicographical comparisons of SWAP gates for some
 * architecture and set of interacting qubits.
 * Used in the 'LexiRoute' method for routing subcircuits as part of the
 * MappingManager framework.
 * Used in solution presented in "On the qubit routing problem" ->
 * arXiv:1902.08091
 */
class LexicographicalComparison {
 public:
  /**
   * Class constructor
   * @param _architecture Architecture object for calcuating distances from
   * @param _interacting_nodes Pairs of physical Node with interacting logical
   * Qubit
   */
  LexicographicalComparison(
      const ArchitecturePtr& _architecture,
      const interacting_nodes_t& _interacting_nodes);

  /**
   * Modifies some distances object by reference.
   * Updates the distance between pair Node in interaction by increment.
   * Increment and Interaction determined by some SWAP.
   *
   * @param distances Distances object updated.
   * @param interaction Node pair increment distance indexing found from
   * @param increment Amount to modify distance index by
   */
  void increment_distances(
      lexicographical_distances_t& distances,
      const std::pair<Node, Node>& interaction, int increment) const;

  /**
   * Returns a held lexicograhically ordered vector of distances between nodes
   * and architectuture class object is constructed from, with changes
   * from increment distances.
   *
   * @return Lexicographically ordered distance vector
   */
  lexicographical_distances_t get_lexicographical_distances() const;

  /**
   * Takes a copy of Distance vector held in object and modifies it to reflect
   * how distance between pairs of interacting nodes in attribute would change
   * given the logical qubits asisgned to the physical node in "swap" swapped.
   *
   * @param swap Physical Node Logical Qubit swapped between to derive copy
   * distance
   */
  lexicographical_distances_t get_updated_distances(const swap_t& swap) const;

  /**
   * For each swap in candidate_swaps, removes swap from set if the distance
   * vector produced by modifying this->lexicographical_distances by said swap
   * is lexicographically smaller to that produced for any other swap. In this
   * way, only swaps with lexicographically identical swap for the given
   * interacting nodes remain after the method is called.
   *
   * @param candidate_swaps Potential pairs of nodes for comparing and removing
   */
  void remove_swaps_lexicographical(swap_set_t& candidate_swaps) const;

 private:
  ArchitecturePtr architecture_;
  lexicographical_distances_t lexicographical_distances;
  interacting_nodes_t interacting_nodes_;
};

}  // namespace tket

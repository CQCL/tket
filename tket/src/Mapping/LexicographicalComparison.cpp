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

#include "Mapping/LexicographicalComparison.hpp"

namespace tket {

LexicographicalComparison::LexicographicalComparison(
    const ArchitecturePtr& _architecture,
    const interacting_nodes_t& _interacting_nodes)
    : architecture_(_architecture), interacting_nodes_(_interacting_nodes) {
  unsigned diameter = this->architecture_->get_diameter();

  lexicographical_distances_t distance_vector(diameter, 0);
  for (const auto& interaction : this->interacting_nodes_) {
    // If Node not in architecture, don't add
    if (!this->architecture_->node_exists(interaction.first) ||
        !this->architecture_->node_exists(interaction.second)) {
      throw LexicographicalComparisonError(
          "Constructor passed some interacting node not in architecture.");
    }
    // key->value already copied, assign reverse to map for later ease
    this->interacting_nodes_[interaction.second] = interaction.first;
    unsigned distance = this->architecture_->get_distance(
        interaction.first, interaction.second);
    if (distance > 0) {
      ++distance_vector[diameter - distance];
    }
  }
  this->lexicographical_distances = distance_vector;
}

void LexicographicalComparison::increment_distances(
    lexicographical_distances_t& distances,
    const std::pair<Node, Node>& interaction, int increment) const {
  const unsigned distances_index =
      this->architecture_->get_diameter() -
      this->architecture_->get_distance(interaction.first, interaction.second);
  if (distances[distances_index] == 0 && increment < 0) {
    throw LexicographicalComparisonError(
        "Negative increment value is larger than value held at index, "
        "modification not allowed.");
  }
  distances[distances_index] += increment;
}

/**
 * getter
 */
lexicographical_distances_t
LexicographicalComparison::get_lexicographical_distances() const {
  return this->lexicographical_distances;
}

/**
 * get_updated_distances
 * updates the "distance vector" (this->lexicographical_distances) to reflect
 * the distance between interacting logical qubits given that the logical qubits
 * present in "swap" have swapped physical qubits (Node)
 */
lexicographical_distances_t LexicographicalComparison::get_updated_distances(
    const swap_t& swap) const {
  // make a copy of base lexicographical distances
  lexicographical_distances_t copy = this->lexicographical_distances;
  if (swap.first == swap.second) {
    return copy;
  }
  auto iq_it = this->interacting_nodes_.find(swap.first);
  // first condition => first node not interacting with self, so update
  // distances
  if (iq_it != this->interacting_nodes_.end()) {
    // update distances due to first swap node and qubit its interating with
    // (assuming swap)
    Node interacting = iq_it->second;
    if (interacting != swap.second) {
      increment_distances(copy, {swap.first, interacting}, -2);
      // updates distances due to second swap node and qubit first is
      // interacting with
      increment_distances(copy, {swap.second, interacting}, 2);
    }
  }
  iq_it = this->interacting_nodes_.find(swap.second);
  // => second node not interacting with self
  if (iq_it != this->interacting_nodes_.end()) {
    Node interacting = iq_it->second;
    if (interacting != swap.first) {
      // update distances due to second node and qubit its interacting with
      increment_distances(copy, {swap.second, interacting}, -2);
      // update distannces due to frist node and qubit second node is
      // interacting with
      increment_distances(copy, {swap.first, interacting}, 2);
    }
  }
  return copy;
}

/**
 * remove_swaps_lexicographical
 * value x at index i of this->lexicographical_distancs => x logical qubits
 * distance (diameter - i) away from the logical qubit they should be
 * interacting with For each swap (swap_t) in "candidate_swaps" a
 * new distances object is created given interacting_qubits Each distance for
 * each swap is lexicographically compared If a distance is lexicographically
 * larger than any other its corresponding swap is removed from candidate_swaps
 * Therefore swaps remaining in candidate_swaps after this process are
 * lexicographically identical for implied logical->physical qubit mapping and
 * interacting logical
 */
void LexicographicalComparison::remove_swaps_lexicographical(
    swap_set_t& candidate_swaps) const {
  auto it = candidate_swaps.begin();
  lexicographical_distances_t winning_distances =
      this->get_updated_distances(*it);
  swap_set_t preserved_swaps = {*it};
  ++it;
  for (; it != candidate_swaps.end(); ++it) {
    lexicographical_distances_t comparison_distances =
        this->get_updated_distances(*it);

    if (comparison_distances < winning_distances) {
      preserved_swaps = {*it};
      winning_distances = comparison_distances;
    } else if (comparison_distances == winning_distances) {
      preserved_swaps.insert(*it);
    }
  }
  candidate_swaps = preserved_swaps;
}
}  // namespace tket
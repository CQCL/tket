// Copyright 2019-2023 Cambridge Quantum Computing
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
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
class RNG;
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** Given an unassigned PV, choose a possible target vertex TV from Dom(PV)
 * so that the search can proceed with trying pv->tv.
 * Different heuristics may give better performance.
 */
class ValueOrdering {
 public:
  ValueOrdering();

  /** Choose TV with larger degrees, with some randomness (e.g.,
   * "Solution-Biased Search").
   * @param possible_values The domain of the pattern vertex PV. This MUST have
   * at least one TV.
   * @param target_ndata Data about graphs, etc. etc. (and also RNG) to assist
   * with the decision.
   * @param rng A random number generator.
   * @return The chosen TV from Dom(PV).
   */
  VertexWSM get_target_value(
      const boost::dynamic_bitset<>& possible_values,
      const NeighboursData& target_ndata, RNG& rng);

 private:
  struct HighDegreeVerticesData {
    std::vector<VertexWSM> vertices;

    // Vertices with higher degree have higher "mass",
    // with the probability of being chosen proportional
    // to the mass.
    unsigned mass;
  };

  // Element [i] is for vertices of degree D-i.
  std::vector<HighDegreeVerticesData> m_entries_for_high_degree_vertices;

  // Fills m_entries_for_high_degree_vertices.
  void fill_entries_for_high_degree_vertices(
      const boost::dynamic_bitset<>& possible_values,
      const NeighboursData& target_ndata);

  // Once fill_entries_for_high_degree_vertices has been called,
  // select a vertex at random,
  // biased so that the probability is proportional to the mass.
  VertexWSM get_random_choice_from_data(RNG& rng) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

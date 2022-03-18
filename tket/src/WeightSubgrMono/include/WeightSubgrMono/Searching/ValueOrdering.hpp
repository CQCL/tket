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
#include <set>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct SharedData;

/** Given an unassigned PV, choose a possible target vertex TV from Dom(PV)
 * so that the search can proceed with trying pv->tv.
 * Different heuristics may give better performance.
 */
class ValueOrdering {
 public:
  ValueOrdering();

  /** Choose TV with larger degrees, with some randomness ("Solution-Biased
   * Search").
   * @param possible_values The domain of the pattern vertex PV.
   * @param shared_data Data about graphs, etc. etc. (and also RNG) to assist
   * with the decision.
   * @return The chosen TV from Dom(PV).
   */
  VertexWSM get_target_value(
      const std::set<VertexWSM>& possible_values, SharedData& shared_data);

 private:
  struct HighDegreeVerticesData {
    std::vector<VertexWSM> vertices;

    // Vertices with higher degree have higher "mass",
    // with the probability of being chosen proportional
    // to the mass.
    unsigned mass;
  };

  // Element [i] is for vertices of degree D-i.
  std::vector<HighDegreeVerticesData> m_data;

  // Fills m_data.
  void fill_data(
      const std::set<VertexWSM>& possible_values, SharedData& shared_data);

  // Once fill_data has been called, select a vertex
  // at random, biased so that the probability is proportional
  // to the mass.
  VertexWSM get_random_choice_from_data(SharedData& shared_data) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

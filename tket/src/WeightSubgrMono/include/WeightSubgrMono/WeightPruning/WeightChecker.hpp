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
#include <memory>
#include <optional>

#include "../GraphTheoretic/GeneralStructs.hpp"
#include "WeightNogoodDetectorManager.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;
class NodeWSM;
class SearchBranch;
class WeightNogoodDetector;

class WeightChecker {
 public:
  /** The weight checking will not take place until the last possible
   * moment, to ensure the best results.
   * When that does happen, it will call the search branch for information
   * (specifically, to get the complete list of used target vertices).
   */
  WeightChecker(
      const NeighboursData& pattern_neighbours_data,
      const NeighboursData& target_neighbours_data,
      const SearchBranch& search_branch, WeightWSM total_p_edge_weights);

  ~WeightChecker();

  struct Result {
    bool nogood;

    /** If not null, we've newly discovered that the target vertex TV
     * is invalid (e.g., there's simply no PV which could ever map to it,
     * without causing the total scalar product to be too high).
     */
    std::optional<VertexWSM> invalid_t_vertex;
  };

  /** Try to find a weight nogood (i.e., an impossible current position).
   * @param possible_assignments Data for all unassigned vertices.
   * @param max_extra_scalar_product The maximum extra weight we allow (the
   * whole point of this nogood detection attempt).
   * @return Information about a possible weight nogood.
   */
  Result operator()(const NodeWSM& node, WeightWSM max_extra_scalar_product);

 private:
  const NeighboursData& m_pattern_neighbours_data;
  const NeighboursData& m_target_neighbours_data;
  const SearchBranch& m_search_branch;
  WeightNogoodDetectorManager m_manager;

  // Will be initialised at the point of first use.
  // Delaying as long as possible can give better results.
  std::unique_ptr<WeightNogoodDetector> m_detector_ptr;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

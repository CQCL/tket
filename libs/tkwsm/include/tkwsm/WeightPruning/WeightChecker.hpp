// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tkwsm/GraphTheoretic/GeneralStructs.hpp"
#include "tkwsm/WeightPruning/WeightNogoodDetectorManager.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;
class NeighboursData;
class SearchBranch;
class WeightNogoodDetector;

/** This wraps both a WeightNogoodDetector object (used to estimate weights,
 * and hence break off early if we can prove that any solution from this point
 * would be rejected, for having too high a scalar product),
 * and a WeightNogoodDetectorManager object
 * (used to monitor progress, and decide whether or not we even want
 * to attempt detection).
 */
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
      const SearchBranch& search_branch, WeightWSM total_p_edge_weights,
      std::set<VertexWSM>& impossible_target_vertices);

  ~WeightChecker();

  /** Try to find a weight nogood (i.e., an impossible current position).
   * @param accessor Object to retrieve data about domains.
   * @param max_extra_scalar_product The maximum extra scalar product allowed,
   * for an acceptable solution (the whole point of this nogood detection
   * attempt).
   * @return False if the current node is a nogood.
   */
  bool check(
      const DomainsAccessor& accessor, WeightWSM max_extra_scalar_product);

  /** Information about many target vertices were/are still valid. */
  struct TVData {
    /** The number of used TV initially passed into the weight detector. */
    std::size_t initial_number_of_tv;

    /** The number of valid TV finally. */
    std::size_t final_number_of_tv;
  };

  /** Null if the weight nogood detector was never activated. */
  std::optional<TVData> get_tv_data_opt();

 private:
  const NeighboursData& m_pattern_neighbours_data;
  const NeighboursData& m_target_neighbours_data;
  const SearchBranch& m_search_branch;
  WeightNogoodDetectorManager m_manager;
  TVData m_tv_data;
  std::set<VertexWSM>& m_impossible_target_vertices;

  // Will be initialised at the point of first use.
  // Delaying as long as possible can give better results.
  std::unique_ptr<WeightNogoodDetector> m_detector_ptr;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

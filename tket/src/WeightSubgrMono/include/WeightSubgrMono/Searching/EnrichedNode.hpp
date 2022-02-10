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
#include <map>
#include <optional>
#include <set>
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"
#include "SearchNodeWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

/** The search node with raw data,
 * together with other associated necessary data for reduction.
 * All reduction objects should autodetect whether
 * further reduction is needed.
 */
struct EnrichedNode {
  SearchNodeWrapper node_wrapper;
  std::size_t n_assignments_processed_by_all_diff_propagator;

  /** All data attached to the SearchNode is cleared,
   * as if it is a new node requiring reduction etc.
   */
  void clear_enriched_data();
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

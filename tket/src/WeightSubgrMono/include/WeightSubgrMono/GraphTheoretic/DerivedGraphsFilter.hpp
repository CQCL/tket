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

#include "DerivedGraphsCalculator.hpp"
#include "DerivedGraphsUpdater.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct DerivedGraphs;
struct FixedData;

/** A vertex filter, i.e. uses properties of derived graphs to check
 * if a single PV->TV assignment may be possible.
 */
class DerivedGraphsFilter {
 public:
  explicit DerivedGraphsFilter(const FixedData& fixed_data);

  ~DerivedGraphsFilter();

  /** Does the individual PV->TV mapping appear to be compatible?
   * @param pv A pattern vertex
   * @param tv A target vertex
   * @param fixed_data Contains necessary data (specifically, NeighboursData),
   * since the information is calculated lazily.
   * @return True if the mapping appears to be possible, false if it is
   * DEFINITELY always impossible, regardless of the other assignments. Thus, if
   * false, TV can be deleted from Dom(PV) at every level.
   */
  bool is_compatible(VertexWSM pv, VertexWSM tv, const FixedData& fixed_data);

  /** For convenience, return the internal derived graphs, to be shared by other
   * objects.
   * @return The derived graphs.
   */
  DerivedGraphs& get_derived_pattern_graphs();
  DerivedGraphs& get_derived_target_graphs();

 private:
  DerivedGraphsCalculator m_calculator;
  DerivedGraphsStorage m_storage;

  std::unique_ptr<DerivedGraphsUpdaterPair> m_updater_pair_ptr;

  // All PV->TV mappings previously calculated to be compatible.
  PossibleAssignments m_compatible_assignments;

  // All PV->TV mappings previously calculated to be impossible.
  PossibleAssignments m_impossible_assignments;

  // For fast checking of PV->TV, we need the edge weights
  // in derived graphs to be sorted (whereas, initially they are stored
  // ordered by vertex instead, for fast neighbour detection).
  // We store them here rather than bothering the graphs objects directly.
  // Note that we only ever need at most 1 object from each map
  // simultaneously, so they are independent and we don't need to worry
  // about invalidated references.
  std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>
      m_d2_pattern_weights;
  std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>
      m_d2_target_weights;

  std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>
      m_d3_pattern_weights;
  std::map<VertexWSM, std::vector<DerivedGraphStructs::Count>>
      m_d3_target_weights;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

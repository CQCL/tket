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
#include "DerivedGraphsContainer.hpp"
#include <forward_list>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;

/** A vertex filter, i.e. uses properties of derived graphs to check
 * if a single PV->TV assignment may be possible.
 */
class DerivedGraphsFilter {
public:

  /** Does the individual PV->TV mapping appear to be compatible?
   * @param pv A pattern vertex
   * @param tv A target vertex
   * @param fixed_data Contains necessary data (specifically, NeighboursData), since the information is calculated lazily.
   * @return True if the mapping appears to be possible, false if it is DEFINITELY always impossible, regardless of the other assignments. Thus, if false, TV can be deleted from Dom(PV) at every level.
   */
  bool is_compatible(VertexWSM pv, VertexWSM tv, const FixedData& fixed_data);

  /** For convenience, return the internal container, to be shared by other objects.
   * @return The internal container with derived graphs data.
  */
  DerivedGraphsContainer& get_container();

private:

  DerivedGraphsContainer m_container;

  // All PV->TV mappings previously calculated to be compatible.
  PossibleAssignments m_compatible_assignments;

  // All PV->TV mappings previously calculated to be impossible.
  PossibleAssignments m_impossible_assignments;
};



}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

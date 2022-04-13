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
#include "../Searching/NodeWSM.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NearNeighboursData;
class NeighboursData;

/** This reduces all domains whenever a SINGLE new pv->tv assignment is made. */
class DistancesRawReducer {
 public:
  DistancesRawReducer(
      const NeighboursData& pattern_ndata,
      NearNeighboursData& pattern_near_ndata,
      const NeighboursData& target_ndata,
      NearNeighboursData& target_near_ndata);

  enum class Result {
    /** The specific PV->TV is impossible, independently of any other values. */
    IMPOSSIBLE_ASSIGNMENT,

    /** After intersecting various domains implied by PV->TV, some domain has
       become empty, so the node is invalid. */
    IMPOSSIBLE_NODE,

    /** The domain intersections for PV->TV have all been applied, and the node
       is not yet invalidated. */
    SUCCESS
  };

  /** Performs the reduction: examines all p-vertices at distance exactly d from
   * pv. */
  // Result operator()(VertexWSM pv, VertexWSM tv, NodeWSM& node, unsigned
  // distance);
  Result operator()(
      const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node,
      unsigned distance);

  /** Checks if pv->tv appears to be valid at first glance. */
  bool check(
      const std::pair<VertexWSM, VertexWSM>& assignment,
      unsigned distance) const;

 private:
  const NeighboursData& m_pattern_ndata;
  NearNeighboursData& m_pattern_near_ndata;
  const NeighboursData& m_target_ndata;
  NearNeighboursData& m_target_near_ndata;

  std::vector<VertexWSM> m_work_vector;
  std::set<VertexWSM> m_work_set;

  /** For distance=1, a special case. */
  // Result reduce_neighbours(VertexWSM pv, VertexWSM tv, NodeWSM& node);
  Result reduce_neighbours(
      const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

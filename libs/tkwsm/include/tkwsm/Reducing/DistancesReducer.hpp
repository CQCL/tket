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
#include "ReducerWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NearNeighboursData;
class NeighboursData;

/** If pv->tv is made, ensures that all pattern vertices v' with
 * distance(pv,v')=d have Domain(v') a subset of
 *    {u : distance(tv,u) <= d}.
 * Note that d >= 2 is required; a separate NeighboursReducer deals
 * with the special case d=1.
 */
class DistancesReducer : public ReducerInterface {
 public:
  DistancesReducer(
      NearNeighboursData& pattern_near_ndata,
      NearNeighboursData& target_near_ndata, unsigned distance);

  virtual bool check(std::pair<VertexWSM, VertexWSM> assignment) override;

  virtual ReductionResult reduce(
      std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
      boost::dynamic_bitset<>& work_set) override;

 private:
  NearNeighboursData& m_pattern_near_ndata;
  NearNeighboursData& m_target_near_ndata;

  // A bit crude to store the distance and have separate objects,
  // but all the heavy data is stored elsewhere and accessed by reference,
  // so actually not inefficient, and it fits the framework better.
  const unsigned m_distance;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

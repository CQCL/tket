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
#include "ReducerWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NeighboursData;

/** This is like DistancesReducer, but with d=1, i.e.
 * simply ensure, when pv->tv is made, that all neighbours of pv
 * have domains contained within the set of neighbours of tv.
 */
class NeighboursReducer : public ReducerInterface {
 public:
  NeighboursReducer(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata);

  virtual bool check(std::pair<VertexWSM, VertexWSM> assignment) override;

  virtual ReductionResult reduce(
      std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
      std::set<VertexWSM>& work_set) override;

 private:
  const NeighboursData& m_pattern_ndata;
  const NeighboursData& m_target_ndata;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

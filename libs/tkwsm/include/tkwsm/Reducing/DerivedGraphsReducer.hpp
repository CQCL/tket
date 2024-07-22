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
#include "tkwsm/GraphTheoretic/DerivedGraphs.hpp"
#include "tkwsm/GraphTheoretic/DerivedGraphsCalculator.hpp"
#include "tkwsm/Reducing/ReducerWrapper.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class DomainsAccessor;

/** Use derived graphs (similar to so-called "supplemental graphs")
 * to restrict domains. There are various graph transformations
 * G -> G' which change edges, but preserve vertices, i.e. V(G)=V(G'),
 * with the following property:
 * whenever a monomorphism  f : V(P) -> V(T)  exists, so by definition f
 * creates a well-defined mapping  F : E(P) -> E(T),  then the new graphs
 * P', T' have the property that the SAME  f : V(P') -> V(T')
 * is a monomorphism from P' to T'.
 *
 * A simple example is the D(n) graph (not standard notation):
 * given any graph G (without edge weights), define G' = D(n)[G] by:
 * (u,v) are adjacent in G'  <==>  there exists a path of length exactly n
 * joining u to v.
 * Furthermore, we can let the edge u--v in D(n) have a weight, equal to
 * the number of distinct paths of length n joining u to v.
 *
 * Note that this is related, but different, to saying  Distance(u,v)=n,
 * which would NOT work. Note also that D(n) would become very expensive
 * to compute for larger n, but we only consider D(2), D(3).
 */
class DerivedGraphsReducer : public ReducerInterface {
 public:
  DerivedGraphsReducer(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata);

  virtual bool check(std::pair<VertexWSM, VertexWSM> assignment) override;

  virtual ReductionResult reduce(
      std::pair<VertexWSM, VertexWSM> assignment, DomainsAccessor& accessor,
      boost::dynamic_bitset<>& work_set) override;

 private:
  // The m_calculator object is shared between the m_derived_pattern_graphs
  // and m_derived_target_graphs objects.
  DerivedGraphsCalculator m_calculator;

  DerivedGraphs m_derived_pattern_graphs;
  DerivedGraphs m_derived_target_graphs;

  // We will call this several times with different derived graphs data.
  ReductionResult reduce_with_derived_data(
      const DerivedGraphStructs::NeighboursAndCounts&
          pattern_derived_neighbours_data,
      const DerivedGraphStructs::NeighboursAndCounts&
          target_derived_neighbours_data,
      VertexWSM root_pattern_vertex, DomainsAccessor& accessor,
      boost::dynamic_bitset<>& work_set);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

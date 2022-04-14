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
#include "../GraphTheoretic/DerivedGraphs.hpp"
#include "../GraphTheoretic/DerivedGraphsCalculator.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class NodeWSM;

class DerivedGraphsReducer {
 public:
  DerivedGraphsReducer(
      const NeighboursData& pattern_ndata, const NeighboursData& target_ndata);

  /** Returns false if the pv->tv assignment is impossible,
   * according to the derived graphs.
   */
  bool check(const std::pair<VertexWSM, VertexWSM>& assignment);

  /** Call at the start when we have a new node, to begin reducing. */
  void clear();

  enum class ReductionResult { SUCCESS, NOGOOD, NEW_ASSIGNMENT };

  /** Automatically remembers which assignments in this node
   * have already been processed; if NEW_ASSIGNMENT is returned,
   * then this will resume from where it left off.
   */
  ReductionResult reduce(NodeWSM& node);

 private:
  DerivedGraphStructs::NeighboursAndCountsStorage m_storage;
  DerivedGraphStructs::SortedCountsStorage m_counts_storage;
  DerivedGraphsCalculator m_calculator;
  DerivedGraphs m_derived_pattern_graphs;
  DerivedGraphs m_derived_target_graphs;

  std::set<VertexWSM> m_new_domain;
  std::size_t m_number_of_assignments_processed;

  // Reduce an individual assignment.
  ReductionResult reduce(
      const std::pair<VertexWSM, VertexWSM>& assignment, NodeWSM& node);

  ReductionResult reduce(
      const DerivedGraphStructs::NeighboursAndCounts& pattern_neighbours,
      const DerivedGraphStructs::NeighboursAndCounts& target_neighbours,
      NodeWSM& node);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

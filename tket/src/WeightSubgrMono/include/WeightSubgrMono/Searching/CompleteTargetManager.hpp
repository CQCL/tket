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
#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
struct SearchNode;

/** Problems with complete target graphs are quite different
 * from "ordinary" target graphs; this is a generalisation of
 * the Travelling Salesman Problem, so e.g. simulated annealing
 * might be better.
 */
class CompleteTargetManager {
 public:
  explicit CompleteTargetManager(const FixedData& fixed_data);

  /** All-in-one variable and valure chooser. Returns the next PV->TV assignment
   * to make. */
  std::pair<VertexWSM, VertexWSM> choose_next_assignment(
      const SearchNode& node, const Assignments& assignments) const;

  /** A replacement for the class VariableOrdering. */
  VertexWSM choose_variable(
      const SearchNode& node, const Assignments& assignments) const;

 private:
  const FixedData& m_fixed_data;

  // Crude: choose PV with LARGEST edge weights sum.
  // (Everything must be assigned, so start with the
  // bad vertex FIRST).
  std::map<VertexWSM, WeightWSM> m_pattern_edge_sums;

  // KEY: target vertex TV
  // VALUE: element[i] is the sum of the smallest i+1 target edge weights
  // joining TV
  std::map<VertexWSM, std::vector<WeightWSM>> m_target_partial_edge_sums;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

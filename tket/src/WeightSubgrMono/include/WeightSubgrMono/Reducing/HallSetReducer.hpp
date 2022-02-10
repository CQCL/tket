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
#include <array>
#include <map>
#include <optional>
#include <set>
#include <string>

#include "../GraphTheoretic/GeneralStructs.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

class SearchBranch;
class SearchNodeWrapper;

/** We have possible mappings [v(i), Dom(v(i))],
 * i.e., for each pattern graph v(i), it must map to one of the target
 * vertices in Dom(v(i)).
 * For each set I of indices, let
 *
 *    U = union_{i in I} Dom(v(i))
 *
 * If there is some I with
 *
 *  | union_{i in I} v(i) | < |U|,
 *
 * then no solution is possible (we have detected a nogood):
 * however we assign, there would be some pair of distinct v(i),v(j)
 * with i,j in I which mapped to the same target vertex.
 * This class tries to detect this, although it does not GUARANTEE
 * to detect this in all cases (it doesn't search through all possible I,
 * because there are exponentially many).
 * Also, if ever we find
 *
 *  | union_{i in I} v(i) | = |U|,
 *
 * Then it's clear that only vertices v(i) for i in I can be assigned
 * to tv in U, and so we may erase every tv in U from every other domain.
 */
class HallSetReducer {
 public:
  /** Try to detect inconsistencies and possible reductions,
   * and carry reductions out if found.
   * @param search_node_wrapper Accesses the internal search node data.
   * @param branch The entire search branch, not just this node; needed just to
   * update "assignments" if reduction occurs.
   * @return False if an inconsistency is detected.
   */
  bool reduce(
      SearchNodeWrapper& search_node_wrapper, SearchBranch& branch) const;

 private:
  // See the paper "A Parallel, Backjumping Subgraph Isomorphism
  // Algorithm using Supplemental Graphs" by Ciaran McCreesh
  // and Patrick Prosser; algorithm 6.
  // However, we'll do a simpler version, just repeatedly reducing
  // until everything stops changing.

  // THe following stored data for reuse
  // is purely to cut down on memory allocation.
  mutable std::set<VertexWSM> m_combined_domains;

  struct VariableData {
    VertexWSM vertex;
    std::size_t domain_size;

    // We'll sort by LARGEST domain size first,
    // so we can pop off the Hall sets as they are discovered.
    bool operator<(const VariableData& other) const;
  };

  mutable std::vector<VariableData> m_domain_sizes_and_vertices;

  // ASSUMING that m_domain_sizes_and_vertices has been filled
  // and sorted, searches for a Hall set at the top
  // (i.e., a set of variables {v(1), ..., v(n)} such that
  //  | union_{1 <= i <= n} Dom(v(i)) | = n).
  bool find_and_remove_top_hall_set_block(
      SearchNodeWrapper& search_node_wrapper, SearchBranch& branch) const;

  // A Hall set was definitely found.
  // Remove it, and resort if necessary
  // (if anything changed).
  bool remove_top_hall_set_block(
      std::size_t hall_set_size, SearchNodeWrapper& search_node_wrapper,
      SearchBranch& branch) const;
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

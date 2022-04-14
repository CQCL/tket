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

/** We have possible mappings [v(i), Dom(v(i))],
 * i.e., for each pattern graph vertex v(i), it must map to one of the target
 * vertices in Dom(v(i)).
 * For each set I of indices, let
 *
 *    U = U(I) = union_{i in I} Dom(v(i)).
 *
 * If there is some I with |I| > |U(I)|,
 * then no solution is possible (we have detected a nogood).
 * However we assign, there would be some pair of distinct v(i), v(j)
 * with i,j in I which mapped to the same target vertex.
 * This class tries to detect this, although it does not GUARANTEE
 * to detect this in all cases (it doesn't search through all possible I,
 * because there are exponentially many).
 * Also, if ever we find |I| = |U(I)|,
 * (when we say that {v(i)} is a "Hall set"),
 * then only vertices v(i) for i in I can be assigned
 * to tv in U, and so we may erase every tv in U from every other domain.
 * Of course this generalises the case |I|=1, which is just an assignment,
 * and thus leads to alldiff propagation. However we only consider
 * |I|>1, and let other routines take care of alldiff propagation.
 */
class HallSetDetector {
 public:
  enum class Status { UNINTERESTING, HALL_SET, NOGOOD };

  /** This will only be of interest if the size of the union is <=
   * the number of PV, so it is a HALL_SET or NOGOOD.
   * The data is needed only when the status is HALL_SET,
   * and it may be inaccurate in other cases.
   */
  struct Result {
    /** If it is a Hall set, this will be sorted to allow binary searches.
     * Otherwise, this can be ignored. */
    std::vector<VertexWSM> pattern_vertices;

    /** This is correct if the status is HALL_SET.
     * Otherwise it can be ignored.
     */
    std::set<VertexWSM> union_of_domains;

    Status status;
  };

  enum class Action { CLEAR_DATA, USE_EXISTING_DATA };

  const Result& get_hall_set(
      const PossibleAssignments& possible_assignments, Action action);

 private:
  // See the paper "A Parallel, Backjumping Subgraph Isomorphism
  // Algorithm using Supplemental Graphs" by Ciaran McCreesh
  // and Patrick Prosser; algorithm 6.
  // The key idea is to sort by domain size, and add vertices
  // with smallest domains first.
  // This makes sense because we're trying to find SMALL unions.

  Result m_result;

  struct Data {
    VertexWSM pv;
    std::size_t domain_size;

    // Will be sorted by domain size in reverse order,
    // i.e. largest size first.
    // Thus we can pop back the discovered Hall sets.
    bool operator<(const Data& other) const;
  };

  std::vector<Data> m_domain_sizes_data;

  // Returns false if a nogood is detected.
  bool fill_new_domain_sizes_data(
      const PossibleAssignments& possible_assignments);

  // Returns false if a nogood is detected.
  bool fill_existing_domain_sizes_data(
      const PossibleAssignments& possible_assignments);

  // m_domain_sizes_data has been filled, sorted, and checked to be nonempty.
  // The m_result vertices have been cleared.
  void fill_result_using_nonempty_domain_sizes_data(
      const PossibleAssignments& possible_assignments);

  // m_domain_sizes_data has been filled.
  // union_of_domains in m_result has already been calculated,
  // and checked to be <= the number of PV used to create it.
  // So the result is definitely interesting. Fill in the PV and status.
  void fill_interesting_result(std::size_t number_of_pv);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

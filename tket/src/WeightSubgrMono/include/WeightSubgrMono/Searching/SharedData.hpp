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
#include <map>
#include <optional>
#include <set>

#include "../Common/RNG.hpp"
#include "../Reducing/CloseVerticesFilter.hpp"
#include "../GraphTheoretic/DerivedGraphsFilter.hpp"
#include "SolutionStorage.hpp"
#include "WeightSubgrMono/Reducing/DerivedGraphsReducer.hpp"

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedData;
class SearchBranch;
class VariableOrdering;
class ValueOrdering;

/** Stores all the data and objects containing algorithms which can be shared
 * between different backtracking searchers, including nogoods, etc. etc.
 * which MAY change as searches progress, conveniently in one place.
 *
 * Also manages search branches (calling move up/move down, etc.)
 */
struct SharedData {
  const FixedData& fixed_data;
  SolutionStorage solution_storage;
  CloseVerticesFilter close_vertices_filter;
  DerivedGraphsFilter derived_graphs_filter;
  DerivedGraphsReducer derived_graphs_reducer;
  RNG rng;

  explicit SharedData(const FixedData& fixed_data);

  /** Initialise a new searcher. The CALLER can create multiple
   * SearchBranch objects, stored elsewhere, and interleave searching
   * (once parallel searching with nogoods is implemented).
   * The CALLER must ensure that branches still exist for the lifetime
   * of this object.
   * @param branch The new search branch, to initialise with (possibly updated)
   * initial data (if parallel searching is implemented, so that previously
   * searched subtrees etc. are removed).
   * @return The result of reducing the initial node with all known information.
   */
  ReductionResult initialise(SearchBranch& branch);

  /** Do one single iteration on the given SearchBranch.
   * Automatically stores nogoods, partial/full solutions, etc. etc.
   * @param branch The search branch to update, with one iteration.
   * @param var_ordering A VariableOrdering object - i.e., choose the next
   * variable (pv) to assign. Different heuristics might give better
   * performance.
   * @param val_ordering A ValueOrdering object - i.e., given pv, choose one of
   * the tv from its domain to assign next. Different heuristics might give
   * better performance.
   * @return FALSE if finished, i.e. no more searching to do anywhere: across
   * ALL search branches, it is complete.
   */
  ReductionResult search(
      SearchBranch& branch, VariableOrdering& var_ordering,
      ValueOrdering& val_ordering);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

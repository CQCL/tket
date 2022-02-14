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

#include "CanonicalRelabelling.hpp"

namespace tket {
namespace tsa_internal {

/** Given a raw vertex->vertex mapping which must be enacted exactly (no empty
 * tokens), attempt to find an optimal or near-optimal result in a table, and
 * handle all vertex back-and-forth relabelling.
 */
class ExactMappingLookup {
 public:
  /** If successful, "swaps" will contain a vector of swaps which performs the
   * desired mapping. */
  struct Result {
    std::vector<Swap> swaps;
    bool success;
    bool too_many_vertices;
  };

  /** The Result object is stored internally. Tries to find a sequence of swaps
   * in the table.
   * @param desired_mapping A (source vertex) -> (target vertex) permutation.
   * @param edges Edges which exist between the vertices (equivalently, the
   * swaps which we are permitted to use). Edges with vertices not appearing in
   * desired_mapping will simply be ignored.
   * @param max_number_of_swaps Stop looking in the table if every possible
   * sequence of swaps in the table which enacts the desired mapping exceeds
   * this length (or doesn't exist at all).
   */
  const Result& operator()(
      const VertexMapping& desired_mapping, const std::vector<Swap>& edges,
      unsigned max_number_of_swaps = 16);

  /** Used for partial mapping lookups; like operator(), but does NOT erase the
   * previous result. Overwrites with a new result if an improvement is found.
   * @param desired_mapping A (source vertex) -> (target vertex) permutation.
   * @param edges Edges which exist between the vertices.
   * @param max_number_of_swaps Stop looking in the table once the swap
   * sequences exceed this length.
   */
  const Result& improve_upon_existing_result(
      const VertexMapping& desired_mapping, const std::vector<Swap>& edges,
      unsigned max_number_of_swaps = 16);

 private:
  Result m_result;
  CanonicalRelabelling m_relabeller;

  /** Attempts to fill m_result, given the relabelling to use.
   * If m_result already has a valid solution (i.e., "success" == true),
   * only fills if the new solution has strictly fewer swaps.
   * @param relabelling_result The result of relabelling, for lookup in the raw
   * table.
   * @param old_edges Edges which exist between the vertices before relabelling.
   * @param max_number_of_swaps Stop looking once the swap sequences exceed this
   * length.
   */
  void fill_result_from_table(
      const CanonicalRelabelling::Result& relabelling_result,
      const std::vector<Swap>& old_edges, unsigned max_number_of_swaps);
};

}  // namespace tsa_internal
}  // namespace tket

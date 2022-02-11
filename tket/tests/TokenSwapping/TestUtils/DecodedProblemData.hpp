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

#include <set>
#include <string>

#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

/** For converting a raw string, representing a fixed sequence of swaps,
 *  into a problem for a TSA. */
struct DecodedProblemData {
  /** This is, of course, a possible SOLUTION to the problem, not part
   *  of the problem input data itself. Since we know at least one solution
   *  (this one), we can compare it with our returned solution
   *  to see how good it is.
   */
  std::vector<Swap> swaps;

  /** The desired source->target mapping for a problem. */
  VertexMapping vertex_mapping;

  size_t number_of_vertices;

  /** Do we require the vertex numbers to be {0,1,2,...,m}, with no gaps? */
  enum class RequireContiguousVertices { YES, NO };

  explicit DecodedProblemData(
      const std::string& str,
      RequireContiguousVertices require_contiguous_vertices =
          RequireContiguousVertices::YES);
};

/** For decoding strings like "1e:2d:3c:4b:5a:69:8:8:9:a:b:c:d:e"
 * as seen in FixedCompleteSolutions, which encode
 * the neighbours of vertices 0,1,2,...,N.
 * Only edges(i,j) with i<j will be returned.
 * @param solution_edges_string A string which encodes the edges of an
 * architecture.
 * @return the set of all edges.
 */
struct DecodedArchitectureData {
  std::set<Swap> edges;

  /** The vertex numbers are contiguous, i.e. 0,1,2,...N for some N. */
  size_t number_of_vertices;

  /** Simply without filling any data. */
  DecodedArchitectureData();

  /** Decodes and fills the data upon construction.
   *  @param solution_edges_string A string which encodes the edges of an
   * architecture.
   */
  explicit DecodedArchitectureData(const std::string& solution_edges_string);
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket

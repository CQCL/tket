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

#include "DecodedProblemData.hpp"
#include "TokenSwapping/VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {
namespace tests {

/** Solves a fixed problem using the current best TSA. */
class BestTsaTester {
 public:
  /** Computes a solution to the problem using our best TSA,
   *  checks it, and returns how many swaps it needed.
   * The edges of the graph are directly taken from the list of swaps in the
   * reference solution.
   * @param data The problem data which was decoded from a string.
   * @return The number of swaps returned by our TSA. The calculated swaps are
   * also checked for correctness.
   */
  size_t get_checked_solution_size(const DecodedProblemData& data);

  /** For problems where the architecture is NOT simply given implicitly
   *  by the swap sequence, so we must also pass in the complete set
   *  of edges, some of which might not appear in the final swaps.
   * @param problem_data The data about a specific problem (calculated swaps,
   * etc.)
   * @param architecture_data Data about the architecture for the problem which
   * is NOT deduced implicitly from the problem data itself (i.e., the edges).
   * @return The number of swaps returned by our TSA. The calculated swaps are
   * also checked for correctness.
   */
  size_t get_checked_solution_size(
      const DecodedProblemData& problem_data,
      const DecodedArchitectureData& architecture_data);

 private:
  SwapList m_raw_swap_list;
  DecodedArchitectureData m_architecture_work_data;
  std::vector<std::pair<unsigned, unsigned>> m_edges_vect;
  VertexMapping m_vertex_mapping_copy;
};

}  // namespace tests
}  // namespace tsa_internal
}  // namespace tket

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

#include "HybridTsa.hpp"
#include "SwapListOptimiser.hpp"
#include "SwapListTableOptimiser.hpp"

namespace tket {

/** This class combines all the different token swapping components together
 * in the best known way to get the best overall end-to-end routine
 * (including different heuristics, parameters etc. whose optimal values
 * are unknown, and require experimentation).
 */
class BestFullTsa : public tsa_internal::PartialTsaInterface {
 public:
  BestFullTsa();

  /** We emphasise that, unlike the general PartialTsaInterface, the solution
   * returned is complete, AND includes all known swap list optimisations.
   * Warning: unlike most PartialTsaInterface objects, the vertex_mapping
   * is NOT updated. (There's no point for a full TSA).
   * @param swaps The list of swaps to append to (does not clear first).
   * @param vertex_mapping The current desired mapping, giving (current source
   * vertex)->(target vertex) mappings. NOT updated at the end.
   * @param distances An object to calculate distances between vertices.
   * @param neighbours An object to calculate adjacent vertices to any given
   * vertex.
   * @param path_finder An object to calculate a shortest path between any
   *    pair of vertices. (Of course, paths might not be unique if the graph
   *    is not a tree).
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours,
      tsa_internal::RiverFlowPathFinder& path_finder) override;

 private:
  tsa_internal::HybridTsa m_hybrid_tsa;
  tsa_internal::SwapListOptimiser m_swap_list_optimiser;
  tsa_internal::SwapListTableOptimiser m_table_optimiser;
};

}  // namespace tket

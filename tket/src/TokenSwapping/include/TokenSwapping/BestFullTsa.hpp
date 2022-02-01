// Copyright 2019-2021 Cambridge Quantum Computing
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

#include "ArchitectureMapping.hpp"
#include "HybridTsa00.hpp"
#include "RNG.hpp"
#include "SwapListOptimiser.hpp"
#include "TableLookup/SwapListTableOptimiser.hpp"

namespace tket {
namespace tsa_internal {

/** To enable easier experimentation, keep this up-to-date with the best
 *  end-to-end known default options, but also make it possible to change
 *  the options.
 *  Also include the best known postprocessing swap list optimisations.
 */
class BestFullTsa : public PartialTsaInterface {
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
      PathFinderInterface& path_finder) override;

  /** Wrapper around the main append_partial_solution function, but constructing
   * and using the best known PathFinderInterface object. The DistancesInterface
   * and NeighboursInterface objects will automatically be constructed.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping. Will be updated with
   * the new added swaps.
   *  @param arch_mapping An ArchitectureMapping object, which knows the graph,
   * how to do Node <-> vertex size_t conversions, etc.
   */
  void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      const ArchitectureMapping& arch_mapping);

  /** For experiments, provide access to the internal stored TSA object. This
   * function may be deleted later!
   *  @return Reference to the internal stored TSA object.
   */
  // HybridTsa00& get_hybrid_tsa_for_testing();

 private:
  HybridTsa00 m_hybrid_tsa;
  SwapListOptimiser m_swap_list_optimiser;
  SwapListTableOptimiser m_table_optimiser;
  RNG m_rng;
};

}  // namespace tsa_internal
}  // namespace tket

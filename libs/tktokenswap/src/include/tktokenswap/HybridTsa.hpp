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

#include "CyclesPartialTsa.hpp"
#include "TrivialTSA.hpp"

namespace tket {
namespace tsa_internal {

/** A full end-to-end TSA, combining the partial cycles TSA
 *  (hopefully good) with the full "trivial" TSA (not so good).
 */
class HybridTsa : public PartialTsaInterface {
 public:
  HybridTsa();

  /** For the current token configuration, calculate a sequence of swaps
   *  to move all tokens home, and append them to the given list.
   *  As this is a full TSA, it guarantees to find a solution.
   *  @param swaps The list of swaps to append to.
   *  @param vertex_mapping The current desired mapping.
   *  @param distances An object to calculate distances between vertices.
   *  @param neighbours An object to calculate adjacent vertices to any given
   * vertex.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices.
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours,
      RiverFlowPathFinder& path_finder) override;

 private:
  CyclesPartialTsa m_cycles_tsa;
  TrivialTSA m_trivial_tsa;
};

}  // namespace tsa_internal
}  // namespace tket

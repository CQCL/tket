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

#include "DistancesInterface.hpp"
#include "NeighboursInterface.hpp"
#include "RiverFlowPathFinder.hpp"
#include "VertexMappingFunctions.hpp"

namespace tket {
namespace tsa_internal {

/** TSA stands for Token Swapping Algorithm.
 *  A "partial TSA" is allowed to give up (not calculate any swaps),
 *  even when the tokens are not all home.
 *  The hope is that different partial TSAs can be combined to give
 *  a good full TSA (i.e., one which always finds a complete solution).
 */
class PartialTsaInterface {
 public:
  /** The algorithm is allowed to fail (not calculate any swaps),
   *  but when it DOES return swaps, it is required to decrease L
   *  (the sum of the distances of each vertex containing a token
   *  from its target vertex).
   *  Thus progress is always nonnegative.
   *  Of course, a complete TSA is a special case.
   *  @param swaps The list of swaps to append to (does not clear first).
   *  @param vertex_mapping The current desired mapping. Each key is the
   *    current vertex where a token is; its value is the target vertex
   *    the token wants to reach. Usually, will be updated upon return to be the
   *    new configuration after performing the swaps.
   *  @param distances An object to calculate distances between vertices.
   *  @param neighbours An object to calculate adjacent vertices to any given
   * vertex.
   *  @param path_finder An object to calculate a shortest path between any
   *    pair of vertices. (Of course, paths might not be unique if the graph
   *    is not a tree, so it is an important part of the heuristics that
   *    the returned paths are fairly "consistent", i.e. "nearby" vertex pairs
   *    should return "nearby" paths).
   */
  virtual void append_partial_solution(
      SwapList& swaps, VertexMapping& vertex_mapping,
      DistancesInterface& distances, NeighboursInterface& neighbours,
      RiverFlowPathFinder& path_finder) = 0;

  /** For debugging purposes, every TSA object has a name.
   *  @return The name of this object (not necessarily unique).
   */
  const std::string& name() const;

 protected:
  std::string m_name;
};

}  // namespace tsa_internal
}  // namespace tket

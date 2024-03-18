// Copyright 2019-2024 Cambridge Quantum Computing
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
#include <tkwsm/GraphTheoretic/GeneralStructs.hpp>

namespace tket {
namespace WeightedSubgraphMonomorphism {

struct FixedArchitectures {
  // Various "heavy hexagon" (or "brick wall pattern") architectures.
  // Sets all weights equal to 1.
  static GraphEdgeWeights get_ibm_brooklyn_65_qubits();
  static GraphEdgeWeights get_ibm_montreal_27_qubits();
  static GraphEdgeWeights get_ibm_guadalupe_16_qubits();
  static GraphEdgeWeights get_ibm_perth_7_qubits();

  // Returns a line with specific vertex labels, and all weights equal to 1.
  static GraphEdgeWeights get_path_qubits(
      const std::vector<VertexWSM>& vertices, bool allow_cycles = true);

  // Returns the path with vertices [i,i+1,i+2,...,j].
  static GraphEdgeWeights get_path_qubits(VertexWSM first, VertexWSM last);

  // Sets all weights equal to 1, and checks that the edges are new.
  static void add_edges(
      GraphEdgeWeights& data, const std::vector<EdgeWSM>& edges,
      bool require_edges_to_be_new = true);

  static void add_paths(
      GraphEdgeWeights& data, const std::vector<std::vector<VertexWSM>>& paths,
      bool allow_cycles = true, bool require_edges_to_be_new = true);

  // Adds the edges from the second to the first, inplace, checking that th
  static void merge_first_with_second(
      GraphEdgeWeights& first, const GraphEdgeWeights& second,
      bool require_edges_to_be_new = true);
};

}  // namespace WeightedSubgraphMonomorphism
}  // namespace tket

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

#include "tket/Placement/Placement.hpp"

namespace tket {

GraphPlacement::Frontier::Frontier(const Circuit& _circ) : circ(_circ) {
  VertexVec input_slice;
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  boolean_in_edges = std::make_shared<b_frontier_t>();

  for (const Qubit& qb : circ.all_qubits()) {
    Vertex input = circ.get_in(qb);
    input_slice.push_back(input);
    Edge candidate = circ.get_nth_out_edge(input, 0);
    quantum_in_edges->insert({qb, circ.skip_irrelevant_edges(candidate)});
  }
  for (const Bit& bit : circ.all_bits()) {
    Vertex input = circ.get_in(bit);
    EdgeVec candidates = circ.get_nth_b_out_bundle(input, 0);
    boolean_in_edges->insert({bit, candidates});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, boolean_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}

void GraphPlacement::Frontier::next_slicefrontier() {
  quantum_in_edges = std::make_shared<unit_frontier_t>();
  boolean_in_edges = std::make_shared<b_frontier_t>();
  for (const std::pair<UnitID, Edge>& pair : quantum_out_edges->get<TagKey>()) {
    Edge new_e = circ.skip_irrelevant_edges(pair.second);
    quantum_in_edges->insert({pair.first, new_e});
    Vertex targ = circ.target(new_e);
    EdgeVec targ_classical_ins =
        circ.get_in_edges_of_type(targ, EdgeType::Boolean);
    boolean_in_edges->insert(
        {Bit("frontier_bit", pair.first.index()), targ_classical_ins});
  }

  CutFrontier next_cut = circ.next_cut(quantum_in_edges, boolean_in_edges);
  slice = next_cut.slice;
  quantum_out_edges = next_cut.u_frontier;
}
}  // namespace tket
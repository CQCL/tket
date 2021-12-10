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

#include "Transform.hpp"

namespace tket {

    // For a gate vertex and an edge, traverse back to find the associated input vertex.
    Vertex get_input_from_vertex_edge(Circuit &circ, const Vertex &current_vertex, const Edge &current_outedge) {
        Edge current_e = current_outedge;
        Vertex current_v = current_vertex; 
        while (true) {
            if (is_initial_q_type(circ.get_OpType_from_Vertex(current_v))) {
                return current_v;
            }
            std::tie(current_v, current_e) = circ.get_prev_pair(current_v, current_e);
        }
    }
    // Assume the circuit only contains CZ gates as the two qubit gate
    Transform Transform::reorder_cz(const ArchitecturePtr& architecture) {
        return Transform([architecture](Circuit &circ) {
            bool success = false;
            BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
                if (circ.get_OpType_from_Vertex(vert) == OpType::CZ)
                {
                    // Check if this operation is valid
                    Vertex input_1 = get_input_from_vertex_edge(circ, vert, circ.get_nth_out_edge(vert, 0));
                    Vertex input_2 = get_input_from_vertex_edge(circ, vert, circ.get_nth_out_edge(vert, 1));
                    UnitID q_1 = circ.get_id_from_in(input_1);
                    UnitID q_2 = circ.get_id_from_in(input_2);
                    if (architecture->valid_operation({Node(q_1), Node(q_2)})) {
                        // Remove the vertex
                        circ.remove_vertex(
                                        vert, Circuit::GraphRewiring::Yes,
                                        Circuit::VertexDeletion::No);
                        // Get input edges for in_q and out_q
                        Edge edge_1 = circ.get_nth_out_edge(input_1, 0);
                        Edge edge_2 = circ.get_nth_out_edge(input_2, 0);
                        // Move the gate to the front
                        circ.rewire(vert, {edge_1, edge_2}, {EdgeType::Quantum, EdgeType::Quantum});
                        success = true;
                    }
                }
            }
            return success;
        });
    }
}   // namespace tket
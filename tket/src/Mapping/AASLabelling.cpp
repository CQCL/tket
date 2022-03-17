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

#include "Mapping/AASLabelling.hpp"

namespace tket {

std::pair<bool, unit_map_t> AASLabellingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  bool found_unplaced_qubit = false;
  bool found_unplaced_ppb = false;

  // search for unplaced qubitto speed up the runtime
  for (Qubit q : mapping_frontier->circuit_.all_qubits()) {
    if (!architecture->node_exists(Node(q))) {
      found_unplaced_qubit = true;
      break;
    }
  }
  if (found_unplaced_qubit) {
    std::shared_ptr<unit_frontier_t> next_frontier =
        frontier_convert_vertport_to_edge(
            mapping_frontier->circuit_, mapping_frontier->linear_boundary);

    CutFrontier next_cut = mapping_frontier->circuit_.next_cut(
        next_frontier, std::make_shared<b_frontier_t>());

    for (const Vertex& v : *next_cut.slice) {
      if (mapping_frontier->circuit_.get_OpType_from_Vertex(v) ==
          OpType::PhasePolyBox) {
        TKET_ASSERT(mapping_frontier->circuit_.is_quantum_node(v));
        Op_ptr op_ptr_ppb =
            mapping_frontier->circuit_.get_Op_ptr_from_Vertex(v);

        for (const Edge& e : mapping_frontier->circuit_.get_in_edges_of_type(
                 v, EdgeType::Quantum)) {
          for (const std::pair<UnitID, Edge>& pair :
               next_frontier->get<TagKey>()) {
            if (pair.second == e) {
              if (!architecture->node_exists(Node(pair.first))) {
                found_unplaced_ppb = true;
              }
            }
          }
        }
      }
    }
  }

  if (!found_unplaced_ppb) {
    return {false, {}};
  } else {
    qubit_vector_t q_vec = mapping_frontier->circuit_.all_qubits();
    unit_map_t qubit_to_nodes_place;
    node_set_t node_set_placed;

    for (Qubit q : q_vec) {
      if (architecture->node_exists(Node(q))) {
        qubit_to_nodes_place.insert({q, Node(q)});
        node_set_placed.insert(Node(q));
      }
    }

    node_vector_t nodes_vec = architecture->get_all_nodes_vec();

    // place all unplaced qubits

    for (Qubit q : q_vec) {
      if (!architecture->node_exists(Node(q))) {
        // found unplaced qubit
        // other checks could be added here to avoid placing unused qubits or
        // qubits that are not in an ppb

        unsigned index_to_use = 0;
        while (node_set_placed.find(nodes_vec[index_to_use]) !=
               node_set_placed.end()) {
          ++index_to_use;
        }
        qubit_to_nodes_place.insert({q, nodes_vec[index_to_use]});
        node_set_placed.insert(nodes_vec[index_to_use]);
        mapping_frontier->update_bimaps(
            mapping_frontier->get_qubit_from_circuit_uid(q),
            nodes_vec[index_to_use]);
      }
    }

    mapping_frontier->update_linear_boundary_uids(qubit_to_nodes_place);
    mapping_frontier->circuit_.rename_units(qubit_to_nodes_place);

    return {true, {}};
  }
}

nlohmann::json AASLabellingMethod::serialize() const {
  nlohmann::json j;
  j["name"] = "AASLabellingMethod";
  return j;
}

AASLabellingMethod AASLabellingMethod::deserialize(
    const nlohmann::json& /*j*/) {
  return AASLabellingMethod();
}

}  // namespace tket

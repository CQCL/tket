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

#include "Mapping/AASRoute.hpp"

namespace tket {

AASRouteRoutingMethod::AASRouteRoutingMethod(
    unsigned _aaslookahead, aas::CNotSynthType _cnotsynthtype)
    : cnotsynthtype_(_cnotsynthtype), aaslookahead_(_aaslookahead) {}

std::pair<bool, unit_map_t> AASRouteRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  std::shared_ptr<unit_frontier_t> next_frontier =
      frontier_convert_vertport_to_edge(
          mapping_frontier->circuit_, mapping_frontier->linear_boundary);

  CutFrontier next_cut = mapping_frontier->circuit_.next_cut(
      next_frontier, std::make_shared<b_frontier_t>());

  // search for ppb in slice
  for (const Vertex& v : *next_cut.slice) {
    if (mapping_frontier->circuit_.get_OpType_from_Vertex(v) ==
        OpType::PhasePolyBox) {
      TKET_ASSERT(mapping_frontier->circuit_.is_quantum_node(v));

      unsigned number_of_qubits = mapping_frontier->circuit_.n_in_edges(v);
      unit_vector_t qubit_vec(number_of_qubits);

      // check all qubits of the ppb if they are placed
      bool box_placed = true;
      for (unsigned i = 0; i < number_of_qubits; ++i) {
        const Edge& e = mapping_frontier->circuit_.get_nth_in_edge(v, i);
        for (const std::pair<UnitID, Edge>& pair :
             next_frontier->get<TagKey>()) {
          if (pair.second == e) {
            qubit_vec[i] = Qubit(pair.first);
            if (!architecture->node_exists(Node(pair.first))) {
              box_placed = false;
            }
          }
        }
      }

      // check that the box we are working on is really placed and the check
      // method has been executed
      // this is imporant if the circuit contains more than one ppb and only one
      // of them is placed

      if (box_placed) {
        // get ppb from op
        Op_ptr op_ptr_ppb =
            mapping_frontier->circuit_.get_Op_ptr_from_Vertex(v);
        const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op_ptr_ppb);

        // Circuit circuit_ppb_place(*ppb.to_circuit());

        // generate aritecture to make sure that the result can be inserted into
        // the given circuit by using  flatten_registers
        auto nodes_vec = architecture->get_all_nodes_vec();
        auto edges_vec = architecture->get_all_edges_vec();

        // create maps from qubits/node to int
        std::map<UnitID, Node> orig_node_to_int_node;
        std::map<UnitID, Node> orig_qubit_to_int_node;

        unsigned id_node = 0;

        unsigned n_nodes = architecture->n_nodes();

        for (Node orig_node : nodes_vec) {
          orig_node_to_int_node.insert({orig_node, Node(n_nodes)});
        }

        for (auto orig_qubit : qubit_vec) {
          orig_node_to_int_node[orig_qubit] = Node(id_node);
          ++id_node;
        }

        for (Node orig_node : nodes_vec) {
          if (orig_node_to_int_node[orig_node] == Node(n_nodes)) {
            orig_node_to_int_node[orig_node] = Node(id_node);
            ++id_node;
          }
        }

        // define new arcitecture with int nodes for ppb
        std::vector<Architecture::Connection> new_con;
        for (auto pair : architecture->get_all_edges_vec()) {
          new_con.push_back(
              {orig_node_to_int_node[UnitID(pair.first)],
               orig_node_to_int_node[UnitID(pair.second)]});
        }

        Architecture new_int_arch = Architecture(new_con);

        TKET_ASSERT(architecture->n_nodes() == new_int_arch.n_nodes());

        Circuit result = aas::phase_poly_synthesis(
            new_int_arch, ppb, aaslookahead_, cnotsynthtype_);

        // make sure the circuit can be inserted
        result.flatten_registers();

        // substitute the ppb vertex in the initial circuit with the routed
        // result
        mapping_frontier->circuit_.substitute(result, v);
        return {true, {}};
      }
    }
  }
  return {false, {}};
}

aas::CNotSynthType AASRouteRoutingMethod::get_cnotsynthtype() const {
  return this->cnotsynthtype_;
}

unsigned AASRouteRoutingMethod::get_aaslookahead() const {
  return this->aaslookahead_;
}

nlohmann::json AASRouteRoutingMethod::serialize() const {
  nlohmann::json j;
  j["aaslookahead"] = this->get_aaslookahead();
  j["cnotsynthtype"] = (unsigned)this->get_cnotsynthtype();
  j["name"] = "AASRouteRoutingMethod";
  return j;
}

AASRouteRoutingMethod AASRouteRoutingMethod::deserialize(
    const nlohmann::json& j) {
  unsigned aaslookahead = j.at("aaslookahead").get<unsigned>();
  aas::CNotSynthType cnotsynthtype =
      (aas::CNotSynthType)j.at("cnotsynthtype").get<unsigned>();
  return AASRouteRoutingMethod(aaslookahead, cnotsynthtype);
}

}  // namespace tket


#include "Mapping/AASRoute.hpp"

namespace tket {

AASRouteRoutingMethod::AASRouteRoutingMethod(
    aas::CNotSynthType cnotsynthtype, unsigned aaslookahead) {
  cnotsynthtype_ = cnotsynthtype;
  aaslookahead_ = aaslookahead;
}

bool AASRouteRoutingMethod::check_method(
    const std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  std::shared_ptr<unit_frontier_t> next_frontier =
      frontier_convert_vertport_to_edge(
          mapping_frontier->circuit_, mapping_frontier->quantum_boundary);

  CutFrontier next_cut = mapping_frontier->circuit_.next_cut(
      next_frontier, std::make_shared<b_frontier_t>());

  for (const Vertex& v : *next_cut.slice) {
    if (mapping_frontier->circuit_.get_OpType_from_Vertex(v) ==
        OpType::PhasePolyBox) {
      Op_ptr op_ptr_ppb = mapping_frontier->circuit_.get_Op_ptr_from_Vertex(v);

      bool box_placed = true;

      for (const Command& com : mapping_frontier->circuit_) {
        Op_ptr op = com.get_op_ptr();
        if (*op_ptr_ppb == *op) {
          unit_vector_t qbs = com.get_args();
          for (auto q : qbs) {
            if (!architecture->node_exists(Node(q))) box_placed = false;
          }
        }
      }

      TKET_ASSERT(mapping_frontier->circuit_.is_quantum_node(v));

      if (box_placed) {
        return true;
      } else {
        // found box which is not placed
      }
    }
  }

  return false;
}

unit_map_t AASRouteRoutingMethod::routing_method(
    std::shared_ptr<MappingFrontier>& mapping_frontier,
    const ArchitecturePtr& architecture) const {
  const Architecture arc(*architecture);

  std::shared_ptr<unit_frontier_t> next_frontier =
      frontier_convert_vertport_to_edge(
          mapping_frontier->circuit_, mapping_frontier->quantum_boundary);

  CutFrontier next_cut = mapping_frontier->circuit_.next_cut(
      next_frontier, std::make_shared<b_frontier_t>());

  for (const Vertex& v : *next_cut.slice) {
    if (mapping_frontier->circuit_.get_OpType_from_Vertex(v) ==
        OpType::PhasePolyBox) {
      TKET_ASSERT(mapping_frontier->circuit_.is_quantum_node(v));
      Op_ptr op_ptr_ppb = mapping_frontier->circuit_.get_Op_ptr_from_Vertex(v);

      unit_vector_t qubit_vec;

      bool box_placed = true;

      for (const Command& com : mapping_frontier->circuit_) {
        Op_ptr op = com.get_op_ptr();
        if (*op_ptr_ppb == *op) {
          qubit_vec = com.get_args();
          for (auto q : qubit_vec) {
            // check if the qubit is placed
            if (!architecture->node_exists(Node(q))) box_placed = false;
          }
        }
      }

      // check that the box we are working on is really placed and the check
      // method has been executed
      // this is imporant if the circuit contains more than one ppb and only one
      // of them is placed

      if (box_placed) {
        const PhasePolyBox& ppb = static_cast<const PhasePolyBox&>(*op_ptr_ppb);

        Circuit circuit_ppb_place(*ppb.to_circuit());

        auto nodes_vec = arc.get_all_nodes_vec();

        auto edges_vec = arc.get_all_edges_vec();

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
        for (auto pair : arc.get_all_edges_vec()) {
          new_con.push_back(
              {orig_node_to_int_node[UnitID(pair.first)],
               orig_node_to_int_node[UnitID(pair.second)]});
        }

        Architecture new_int_arch = Architecture(new_con);

        TKET_ASSERT(arc.n_nodes() == new_int_arch.n_nodes());

        Circuit result = aas::phase_poly_synthesis(
            new_int_arch, ppb, aaslookahead_, cnotsynthtype_);

        // make sure the circuit can be inserted
        result.flatten_registers();

        // substitute the ppb vertex in the initial circuit with the routed
        // result
        mapping_frontier->circuit_.substitute(result, v);

      } else {
        TKET_ASSERT(!"box is not placed\n");
      }
    }
  }

  return {};
}

}  // namespace tket

#include "Circuit/CircPool.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

enum class PortType { In, Out };
typedef std::pair<VertPort, PortType> TypedVertPort;
typedef boost::bimap<zx::ZXVert, Vertex> BoundaryVertMap;

bool is_spiderless_optype(const OpType& optype) {
  return optype == OpType::Barrier || optype == OpType::SWAP ||
         optype == OpType::noop;
}

BoundaryVertMap circuit_to_zx_recursive(
    const Circuit& circ, zx::ZXDiagram& zxd, bool add_boundary) {
  // TODO: how can users know the boundary mapping?
  sequenced_map_t<TypedVertPort, zx::ZXVert> vert_lookup;
  BoundaryVertMap bmap;

  // Convert each vertex to ZXDiagram, raise error if not supported
  BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
    // We currently throw an error if the vertex is either
    // 1. conditional, classical, flow
    Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
    if (is_flowop_type(op->get_type()) || is_classical_type(op->get_type()) ||
        op->get_type() == OpType::Conditional) {
      throw Unsupported(
          "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
    switch (op->get_type()) {
      case OpType::Input: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Input, zx::QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::Output: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Output, zx::QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClInput: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Input, zx::QuantumType::Classical);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClOutput: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Output, zx::QuantumType::Classical);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      // Spiderless ops are handled during vertex wiring
      case OpType::Barrier:
      case OpType::noop:
      case OpType::SWAP: {
        continue;
      }
      case OpType::H: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Hbox, zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        zxd.multiply_scalar(0.5);
        break;
      }
      case OpType::Rz: {
        zx::ZXVert zx_vert = zxd.add_vertex(
            zx::ZXType::ZSpider, op->get_params()[0], zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Rx: {
        zx::ZXVert zx_vert = zxd.add_vertex(
            zx::ZXType::XSpider, op->get_params()[0], zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::X: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::XSpider, 1, zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Z: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 1, zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::CX: {
        zx::ZXVert zx_x_vert =
            zxd.add_vertex(zx::ZXType::XSpider, 0, zx::QuantumType::Quantum);
        zx::ZXVert zx_z_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Quantum);
        zxd.add_wire(zx_x_vert, zx_z_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_z_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_z_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_x_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_x_vert});
        zxd.multiply_scalar(2);
        break;
      }
      case OpType::CZ: {
        zx::ZXVert zx_za_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Quantum);
        zx::ZXVert zx_zb_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Quantum);
        zxd.add_wire(zx_za_vert, zx_zb_vert, zx::ZXWireType::H);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_za_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_za_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_zb_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_zb_vert});
        break;
      }
      case OpType::Measure: {
        // Add a decoherence node
        zx::ZXVert zx_measure_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
        // Add a delete operator
        zx::ZXVert zx_delete_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_measure_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_measure_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_delete_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_measure_vert});
        break;
      }
      case OpType::Reset: {
        // Discard
        zx::ZXVert zx_discard_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
        // Add a node to prepare |0>
        zx::ZXVert zx_reset_vert =
            zxd.add_vertex(zx::ZXType::XSpider, 0, zx::QuantumType::Quantum);
        zxd.multiply_scalar(0.5);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_discard_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_reset_vert});
        break;
      }
      case OpType::Collapse: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, zx::QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Create: {
        zx::ZXVert zx_init_vert =
            zxd.add_vertex(zx::ZXType::XSpider, 0, zx::QuantumType::Quantum);
        zxd.multiply_scalar(0.5);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_init_vert});
        break;
      }
      case OpType::Discard: {
        zx::ZXVert zx_discard_vert =
            zxd.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_discard_vert});
        break;
      }
      default:
        // TODO Rebase quantum gates
        if (op->get_type() == OpType::Conditional) {
          // Add a multi-leg swap
          // Add a op or a box
          // Add a multi-leg swap
        }
        if (is_box_type(op->get_type())) {
          const Box& b = static_cast<const Box&>(*op);
          Circuit replacement = *b.to_circuit();
          BoundaryVertMap box_bm =
              circuit_to_zx_recursive(replacement, zxd, false);
          // Find the boundaries of the replacement circuit associated with each
          // port
          EdgeVec q_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Quantum);
          EdgeVec q_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Quantum);
          EdgeVec c_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Classical);
          EdgeVec c_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Classical);
          // TODO:: handle boolean wires
          // std::cout<<"\nA\n";
          for (unsigned i = 0; i < q_in_holes.size(); i++) {
            port_t port = circ.get_target_port(q_in_holes[i]);
            Vertex inp = replacement.get_in(Qubit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::In}, box_bm.right.find(inp)->second});
          }
          for (unsigned i = 0; i < q_out_holes.size(); i++) {
            port_t port = circ.get_source_port(q_out_holes[i]);
            Vertex outp = replacement.get_out(Qubit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::Out},
                 box_bm.right.find(outp)->second});
          }
          for (unsigned i = 0; i < c_in_holes.size(); i++) {
            port_t port = circ.get_target_port(c_in_holes[i]);
            Vertex inp = replacement.get_in(Bit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::In}, box_bm.right.find(inp)->second});
          }
          for (unsigned i = 0; i < c_out_holes.size(); i++) {
            port_t port = circ.get_source_port(c_out_holes[i]);
            Vertex outp = replacement.get_out(Bit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::Out},
                 box_bm.right.find(outp)->second});
          }
        } else {
          throw Unsupported(
              "Cannot convert gate type: " + op->get_name() +
              " to a ZX node, try rebase the gates to use Rx, Rz, X, Z, H, CZ "
              "or CX. \n");
        }
    }
  }

  // Connect the ZX nodes
  BGL_FORALL_EDGES(edge, circ.dag, DAG) {
    Vertex v_s = circ.source(edge);
    Vertex v_t = circ.target(edge);
    port_t p_s = circ.get_source_port(edge);
    port_t p_t = circ.get_target_port(edge);
    // Handle Spiderless ops
    if (is_spiderless_optype(circ.get_OpType_from_Vertex(v_s))) {
      // We only handle in-edges
      continue;
    }
    if (is_spiderless_optype(circ.get_OpType_from_Vertex(v_t))) {
      // Traverse the path to find the next non-spiderless op
      Edge next_e = edge;
      do {
        if (circ.get_OpType_from_Vertex(v_t) == OpType::SWAP) {
          next_e = circ.get_nth_out_edge(
              v_t, (circ.get_target_port(next_e) + 1) % 2);
          v_t = circ.target(next_e);
        } else {
          std::tie(v_t, next_e) = circ.get_next_pair(v_t, next_e);
        }
      } while (is_spiderless_optype(circ.get_OpType_from_Vertex(v_t)));
      p_t = circ.get_target_port(next_e);
    }

    auto it_s = vert_lookup.get<TagKey>().find(
        TypedVertPort(VertPort(v_s, p_s), PortType::Out));
    auto it_t = vert_lookup.get<TagKey>().find(
        TypedVertPort(VertPort(v_t, p_t), PortType::In));
    if (circ.get_edgetype(edge) == EdgeType::Quantum) {
      zxd.add_wire(it_s->second, it_t->second);
    } else {
      zxd.add_wire(
          it_s->second, it_t->second, zx::ZXWireType::Basic,
          zx::QuantumType::Classical);
    }
  }
  return bmap;
}

zx::ZXDiagram circuit_to_zx(const Circuit& circ) {
  zx::ZXDiagram zxd;
  BoundaryVertMap bmap = circuit_to_zx_recursive(circ, zxd, true);
  // TODO return the map as well
  zx::ZXVertVec true_boundary = zxd.get_boundary();
  zx::ZXVertIterator vi, vi_end, next;
  tie(vi, vi_end) = boost::vertices(*zxd.get_graph());
  for (next = vi; vi != vi_end; vi = next) {
    ++next;
    if ((zxd.get_zxtype(*vi) == zx::ZXType::Input ||
         zxd.get_zxtype(*vi) == zx::ZXType::Output) &&
        std::find(true_boundary.begin(), true_boundary.end(), *vi) ==
            true_boundary.end()) {
      zx::WireVec adj_wires = zxd.adj_wires(*vi);
      TKET_ASSERT(adj_wires.size() == 2);
      TKET_ASSERT(zxd.get_qtype(adj_wires[0]) == zxd.get_qtype(adj_wires[1]));
      TKET_ASSERT(
          zxd.get_wire_type(adj_wires[0]) == zxd.get_wire_type(adj_wires[1]));
      zx::ZXVertVec neighbours = zxd.neighbours(*vi);
      zxd.add_wire(
          neighbours[0], neighbours[1], zxd.get_wire_type(adj_wires[0]),
          zxd.get_qtype(adj_wires[0]));
      zxd.remove_vertex(*vi);
    }
  }
  return zxd;
}

}  // namespace tket
#include "Circuit/CircPool.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

using namespace zx;

enum class PortType { In, Out };
typedef std::pair<VertPort, PortType> TypedVertPort;
typedef boost::bimap<ZXVert, Vertex> BoundaryVertMap;

bool is_spiderless_optype(const OpType& optype) {
  return optype == OpType::Barrier || optype == OpType::SWAP ||
         optype == OpType::noop;
}

// A swicth with the on-state controlled by the on_value
std::pair<ZXVert, ZXVert> add_switch(
    ZXDiagram& zxd, const bool& on_value, const QuantumType& qtype) {
  // TODO: what is the scalar for a classical switch?
  zxd.multiply_scalar(std::sqrt(2.));
  ZXVert triangle = zxd.add_vertex(ZXType::Triangle, qtype);
  ZXVert x = zxd.add_vertex(ZXType::XSpider, 0, qtype);
  zxd.add_wire(triangle, x, ZXWireType::Basic, qtype, 1);
  if (on_value) {
    // Turn on the switch if the input is 1
    ZXVert negate = zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Classical);
    zxd.add_wire(triangle, negate, ZXWireType::Basic, qtype, 0);
    return {negate, x};
  }
  return {triangle, x};
}

// --o--s--o--
//   |     |
//   n     n
//   |--G--|
// s, n are switches, G is the conditional zx diagram specified by 'left' and
// 'right'. s is on when the input is 0, n is on when the input is 1
std::pair<std::pair<ZXVert, ZXVert>, std::vector<ZXVert>> add_conditional_zx(
    ZXDiagram& zxd, const ZXVert& left, const ZXVert& right,
    const QuantumType& qtype) {
  ZXVert in = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert out = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  auto switch_i = add_switch(zxd, false, qtype);
  auto switch_s0 = add_switch(zxd, true, qtype);
  auto switch_s1 = add_switch(zxd, true, qtype);
  zxd.add_wire(in, switch_i.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out, switch_i.second, ZXWireType::Basic, qtype);
  zxd.add_wire(in, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(left, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(right, switch_s1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out, switch_s1.second, ZXWireType::Basic, qtype);
  return {{in, out}, {switch_i.first, switch_s0.first, switch_s1.first}};
}

BoundaryVertMap circuit_to_zx_recursive(
    const Circuit& circ, ZXDiagram& zxd, bool add_boundary) {
  // TODO: how can users know the boundary mapping?
  std::map<TypedVertPort, ZXVert> vert_lookup;
  std::map<VertPort, ZXVert> boolean_outport_lookup;
  BoundaryVertMap bmap;

  // Convert each vertex to ZXDiagram, raise error if not supported
  BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
    // We currently throw an error if the vertex is either
    // 1. conditional, classical, flow
    Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
    if (is_flowop_type(op->get_type()) || is_classical_type(op->get_type())) {
      throw Unsupported(
          "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
    switch (op->get_type()) {
      case OpType::Input: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Input, QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::Output: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Output, QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClInput: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Input, QuantumType::Classical);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClOutput: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Output, QuantumType::Classical);
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
        ZXVert zx_vert = zxd.add_vertex(ZXType::Hbox, QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        zxd.multiply_scalar(0.5);
        break;
      }
      case OpType::Rz: {
        ZXVert zx_vert = zxd.add_vertex(
            ZXType::ZSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Rx: {
        ZXVert zx_vert = zxd.add_vertex(
            ZXType::XSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::X: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Z: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::ZSpider, 1, QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::CX: {
        ZXVert zx_x_vert =
            zxd.add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        ZXVert zx_z_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        zxd.add_wire(zx_x_vert, zx_z_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_z_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_z_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_x_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_x_vert});
        zxd.multiply_scalar(2);
        break;
      }
      case OpType::CZ: {
        ZXVert zx_za_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        ZXVert zx_zb_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        zxd.add_wire(zx_za_vert, zx_zb_vert, ZXWireType::H);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_za_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_za_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_zb_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_zb_vert});
        break;
      }
      case OpType::Measure: {
        // Add a decoherence node
        ZXVert zx_measure_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        // Add a delete operator
        ZXVert zx_delete_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_measure_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_measure_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_delete_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_measure_vert});
        break;
      }
      case OpType::Reset: {
        // Discard
        ZXVert zx_discard_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        // Add a node to prepare |0>
        ZXVert zx_reset_vert =
            zxd.add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        zxd.multiply_scalar(0.5);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_discard_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_reset_vert});
        break;
      }
      case OpType::Collapse: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::ZSpider, QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Create: {
        ZXVert zx_init_vert =
            zxd.add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        zxd.multiply_scalar(0.5);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_init_vert});
        break;
      }
      case OpType::Discard: {
        ZXVert zx_discard_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_discard_vert});
        break;
      }
      default:
        // TODO Rebase quantum gates
        if (op->get_type() == OpType::Conditional) {
          const Conditional& cond = static_cast<const Conditional&>(*op);
          Op_ptr inner_op = cond.get_op();
          // First, convert the underlying op to zx.
          op_signature_t inner_sig = inner_op->get_signature();
          Circuit replacement;
          if (is_box_type(inner_op->get_type())) {
            const Box& b = static_cast<const Box&>(*inner_op);
            replacement = *b.to_circuit();
          } else {
            // Assume the inner op doesn't have boolean edges
            unit_vector_t args;
            unsigned q_index = 0;
            unsigned c_index = 0;
            for (const EdgeType& etype : inner_sig) {
              if (etype == EdgeType::Quantum) {
                args.push_back(Qubit(q_index++));
              } else if (etype == EdgeType::Classical) {
                args.push_back(Bit(c_index++));
              }
            }
            replacement = Circuit(q_index, c_index);
            replacement.add_op(inner_op, args);
          }
          BoundaryVertMap box_bm =
              circuit_to_zx_recursive(replacement, zxd, false);
          EdgeVec b_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Boolean);
          // Add a control spider
          ZXVert zx_control_vert =
              zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
          // For each qubit path in the inner op, convert it into a conditional
          // zx.
          unsigned q_index = 0;
          unsigned c_index = 0;
          for (unsigned i = 0; i < inner_sig.size(); i++) {
            QuantumType qtype;
            UnitID reg;
            if (inner_sig[i] == EdgeType::Quantum) {
              qtype = QuantumType::Quantum;
              reg = Qubit(q_index++);
            } else {
              qtype = QuantumType::Classical;
              reg = Bit(c_index++);
            }
            std::pair<ZXVert, ZXVert> boundary;
            ZXVertVec controls;
            Vertex inp = replacement.get_in(reg);
            ZXVert zx_inp = box_bm.right.find(inp)->second;
            Vertex outp = replacement.get_out(reg);
            ZXVert zx_outp = box_bm.right.find(outp)->second;
            std::tie(boundary, controls) =
                add_conditional_zx(zxd, zx_inp, zx_outp, qtype);
            for (const ZXVert& ctr : controls) {
              if (zxd.get_zxtype(ctr) == ZXType::Triangle) {
                zxd.add_wire(
                    ctr, zx_control_vert, ZXWireType::Basic,
                    QuantumType::Quantum, 0);
              }
              else {
                zxd.add_wire(
                    ctr, zx_control_vert, ZXWireType::Basic,
                    QuantumType::Classical);
              }

            }
            // Update lookup
            vert_lookup.insert(
                {{{vert, i + cond.get_width()}, PortType::In}, boundary.first});
            vert_lookup.insert(
                {{{vert, i + cond.get_width()}, PortType::Out},
                 boundary.second});
          }
          // Deal with the conditional bits
          unsigned value = cond.get_value();
          for (unsigned i = 0; i < cond.get_width(); i++) {
            port_t p_t = circ.get_target_port(b_in_holes[i]);
            bool set = value & (1 << (cond.get_width() - i - 1));
            if (!set) {
              ZXVert x_vert =
                  zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Classical);
              zxd.add_wire(
                  x_vert, zx_control_vert, ZXWireType::Basic,
                  QuantumType::Classical);
              vert_lookup.insert({{{vert, p_t}, PortType::In}, x_vert});
            } else {
              vert_lookup.insert(
                  {{{vert, p_t}, PortType::In}, zx_control_vert});
            }
            // Each boolean edge shares a source port with other
            // classical/boolean edges Use the classical z spider to explicitly
            // implement this copy operation
            ZXVert copy_vert =
                zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
            port_t p_s = circ.get_source_port(b_in_holes[i]);
            Vertex vert_s = circ.source(b_in_holes[i]);
            // During the wiring stage, each edge originated from {vert_s, p_s}
            // needs go through copy_vert
            boolean_outport_lookup.insert({{vert_s, p_s}, copy_vert});
          }
        } else if (is_box_type(op->get_type())) {
          const Box& b = static_cast<const Box&>(*op);
          Circuit replacement = *b.to_circuit();
          // First, add the converted box to diagram.
          BoundaryVertMap box_bm =
              circuit_to_zx_recursive(replacement, zxd, false);
          // Then, map the vertports in the box boundary to zx nodes.
          // Assume that a box can't have Boolean input edges, and all Boolean
          // output edges share ports with Classical edges. Therefore we don't
          // have to map Boolean vertports.
          EdgeVec q_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Quantum);
          EdgeVec q_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Quantum);
          EdgeVec c_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Classical);
          EdgeVec c_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Classical);

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
              " to a ZX node, try rebase the gates to use Rx, Rz, "
              "X, Z, H, CZ "
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

    auto it_s =
        vert_lookup.find(TypedVertPort(VertPort(v_s, p_s), PortType::Out));
    auto it_t =
        vert_lookup.find(TypedVertPort(VertPort(v_t, p_t), PortType::In));
    if (circ.get_edgetype(edge) == EdgeType::Quantum) {
      zxd.add_wire(it_s->second, it_t->second);
    } else {
      auto bool_it = boolean_outport_lookup.find(VertPort(v_s, p_s));
      if (bool_it == boolean_outport_lookup.end()) {
        zxd.add_wire(
            it_s->second, it_t->second, ZXWireType::Basic,
            QuantumType::Classical);
      } else {
        if (zxd.degree(bool_it->second) == 0) {
          zxd.add_wire(
              it_s->second, bool_it->second, ZXWireType::Basic,
              QuantumType::Classical);
        }
        zxd.add_wire(
            bool_it->second, it_t->second, ZXWireType::Basic,
            QuantumType::Classical);
      }
    }
  }
  return bmap;
}

ZXDiagram circuit_to_zx(const Circuit& circ) {
  ZXDiagram zxd;
  BoundaryVertMap bmap = circuit_to_zx_recursive(circ, zxd, true);
  // TODO return the map as well
  ZXVertVec true_boundary = zxd.get_boundary();
  ZXVertIterator vi, vi_end, next;
  tie(vi, vi_end) = boost::vertices(*zxd.get_graph());
  for (next = vi; vi != vi_end; vi = next) {
    ++next;
    if ((zxd.get_zxtype(*vi) == ZXType::Input ||
         zxd.get_zxtype(*vi) == ZXType::Output) &&
        std::find(true_boundary.begin(), true_boundary.end(), *vi) ==
            true_boundary.end()) {
      WireVec adj_wires = zxd.adj_wires(*vi);
      TKET_ASSERT(adj_wires.size() == 2);
      TKET_ASSERT(zxd.get_qtype(adj_wires[0]) == zxd.get_qtype(adj_wires[1]));
      TKET_ASSERT(
          zxd.get_wire_type(adj_wires[0]) == zxd.get_wire_type(adj_wires[1]));
      ZXVertVec neighbours = zxd.neighbours(*vi);
      zxd.add_wire(
          neighbours[0], neighbours[1], zxd.get_wire_type(adj_wires[0]),
          zxd.get_qtype(adj_wires[0]));
      zxd.remove_vertex(*vi);
    }
  }
  return zxd;
}

}  // namespace tket
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
  zxd.add_wire(triangle, x, ZXWireType::Basic, qtype);
  if (on_value) {
    // Turn on the switch if the input is 1
    ZXVert negate = zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Classical);
    zxd.add_wire(triangle, negate, ZXWireType::Basic, qtype);
    return {negate, x};
  }
  return {triangle, x};
}

std::pair<ZXVertVec, ZXVertVec> add_conditional_swap(
    ZXDiagram& zxd, const QuantumType& qtype) {
  ZXVert in0 = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert in1 = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert out0 = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert out1 = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  auto switch_i0 = add_switch(zxd, false, qtype);
  auto switch_i1 = add_switch(zxd, false, qtype);
  auto switch_s0 = add_switch(zxd, true, qtype);
  auto switch_s1 = add_switch(zxd, true, qtype);
  zxd.add_wire(in0, switch_i0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out0, switch_i0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(in1, switch_i1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out1, switch_i1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(in0, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out1, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(in1, switch_s1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out0, switch_s1.second, ZXWireType::Basic, qtype);
  return {
      {in0, in1, out0, out1},
      {switch_i0.first, switch_i1.first, switch_s0.first, switch_s1.first}};
}

BoundaryVertMap circuit_to_zx_recursive(
    const Circuit& circ, ZXDiagram& zxd, bool add_boundary) {
  // TODO: how can users know the boundary mapping?
  sequenced_map_t<TypedVertPort, ZXVert> vert_lookup;
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
          EdgeVec ins =
              circ.get_in_edges(vert);
          EdgeVec outs =
              circ.get_all_out_edges(vert);
          EdgeVec q_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Quantum);
          EdgeVec c_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Classical);
          EdgeVec b_in_holes =
              circ.get_in_edges_of_type(vert, EdgeType::Boolean);
          EdgeVec b_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Boolean);
          EdgeVec c_out_holes =
              circ.get_out_edges_of_type(vert, EdgeType::Classical);
          // Add a control spider
          ZXVert zx_control_vert =
              zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);

          for (unsigned i = 0; i < q_in_holes.size(); i++) {
            ZXVertVec swap_boundary0, swap_boundary1, swap_ctrs0, swap_ctrs1;
            // Add the first swap
            std::tie(swap_boundary0, swap_ctrs0) =
                add_conditional_swap(zxd, QuantumType::Quantum);
            Vertex inp = replacement.get_in(Qubit(i));
            ZXVert zx_inp = box_bm.right.find(inp)->second;
            // Connect the bottom wire to the op
            zxd.add_wire(swap_boundary0[3], zx_inp);
            // Add the second swap
            std::tie(swap_boundary1, swap_ctrs1) =
                add_conditional_swap(zxd, QuantumType::Quantum);
            Vertex outp = replacement.get_out(Qubit(i));
            ZXVert zx_outp = box_bm.right.find(outp)->second;
            // Connect the bottom wire to the op
            zxd.add_wire(swap_boundary1[1], zx_outp);
            // Connect the top wires
            zxd.add_wire(swap_boundary0[2], swap_boundary1[0]);
            // Connect switches to control
            for (const ZXVert& sv : swap_ctrs0) {
              zxd.add_wire(sv, zx_control_vert);
            }
            for (const ZXVert& sv : swap_ctrs1) {
              zxd.add_wire(sv, zx_control_vert);
            }
            // Update lookup
            port_t port = circ.get_target_port(q_in_holes[i]);
            vert_lookup.insert(
                {{{vert, port}, PortType::In}, swap_boundary0[0]});
            vert_lookup.insert(
                {{{vert, port}, PortType::Out}, swap_boundary1[2]});
          }
          std::cout<<"\nE\n";
          for (unsigned i = 0; i < c_in_holes.size(); i++) {
            unsigned j = i;
            ZXVertVec swap_boundary0, swap_boundary1, swap_ctrs0, swap_ctrs1;
            // Add the first swap
            std::tie(swap_boundary0, swap_ctrs0) =
                add_conditional_swap(zxd, QuantumType::Classical);
            Vertex inp = replacement.get_in(Bit(j));
            ZXVert zx_inp = box_bm.right.find(inp)->second;
            // Connect the bottom wire to the op
            zxd.add_wire(
                swap_boundary0[3], zx_inp, ZXWireType::Basic,
                QuantumType::Classical);
            // Add the second swap
            std::tie(swap_boundary1, swap_ctrs1) =
                add_conditional_swap(zxd, QuantumType::Classical);
            Vertex outp = replacement.get_out(Bit(j));
            ZXVert zx_outp = box_bm.right.find(outp)->second;
            // Connect the bottom wire to the op
            zxd.add_wire(
                swap_boundary1[1], zx_outp, ZXWireType::Basic,
                QuantumType::Classical);
            // Connect the top wires
            zxd.add_wire(
                swap_boundary0[2], swap_boundary1[0], ZXWireType::Basic,
                QuantumType::Classical);
            // Connect switches to control
            for (const ZXVert& sv : swap_ctrs0) {
              zxd.add_wire(sv, zx_control_vert);
            }
            for (const ZXVert& sv : swap_ctrs1) {
              zxd.add_wire(sv, zx_control_vert);
            }
            // Update lookup
            port_t port = circ.get_target_port(c_in_holes[i]);
            vert_lookup.insert(
                {{{vert, port}, PortType::In}, swap_boundary0[0]});
            vert_lookup.insert(
                {{{vert, port}, PortType::Out}, swap_boundary1[2]});
          }
          std::cout<<"\nF\n";
          // Deal with the conditional bits
          unsigned value = cond.get_value();
          for (unsigned i = 0; i < cond.get_width(); i++) {
            ZXVert z_vert =
                zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
            bool set = value & (1 << (cond.get_width() - i - 1));
            std::cout<<"\nSet:"<<set<<"\n";
            if (!set) {
              ZXVert x_vert =
                  zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Classical);
              zxd.add_wire(
                  z_vert, x_vert, ZXWireType::Basic, QuantumType::Classical);
              zxd.add_wire(
                  x_vert, zx_control_vert, ZXWireType::Basic,
                  QuantumType::Classical);
            } else {
              zxd.add_wire(
                  z_vert, zx_control_vert, ZXWireType::Basic,
                  QuantumType::Classical);
            }
            std::cout<<"\nH1\n"<<b_in_holes.size();
            port_t port = circ.get_target_port(b_in_holes[i]);
            std::cout<<"\nH2\n";
            vert_lookup.insert({{{vert, port}, PortType::In}, z_vert});
          }
          std::cout<<"\nG\n";
        }
        else if (is_box_type(op->get_type())) {
          const Box& b = static_cast<const Box&>(*op);
          Circuit replacement = *b.to_circuit();
          // First, add the converted box to diagram.
          BoundaryVertMap box_bm =
              circuit_to_zx_recursive(replacement, zxd, false);
          // Then, map the vertports in the box boundary to zx nodes.
          // Assume that a box can't have Boolean input edges, and all Boolean output edges
          // share ports with Classical edges. Therefore we don't have to map Boolean vertports.
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

    auto it_s = vert_lookup.get<TagKey>().find(
        TypedVertPort(VertPort(v_s, p_s), PortType::Out));
    auto it_t = vert_lookup.get<TagKey>().find(
        TypedVertPort(VertPort(v_t, p_t), PortType::In));
    if (circ.get_edgetype(edge) == EdgeType::Quantum) {
      zxd.add_wire(it_s->second, it_t->second);
    } else {
      zxd.add_wire(
          it_s->second, it_t->second, ZXWireType::Basic,
          QuantumType::Classical);
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
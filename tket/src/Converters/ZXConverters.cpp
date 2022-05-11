#include "Circuit/CircPool.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

using namespace zx;

enum class PortType { In, Out };
typedef std::pair<VertPort, PortType> TypedVertPort;
typedef std::pair<ZXVert, std::optional<unsigned>> ZXVertPort;
typedef std::vector<ZXVertPort> ZXVertPortVec;
typedef boost::bimap<ZXVert, Vertex> BoundaryVertMap;

bool is_spiderless_optype(const OpType& optype) {
  return optype == OpType::Barrier || optype == OpType::SWAP ||
         optype == OpType::noop;
}

// Add a swicth with the on-state controlled by the on_value.
// qtype indicates the quantum type of the switch.
std::pair<ZXVert, ZXVert> add_switch(
    ZXDiagram& zxd, const bool& on_value, const QuantumType& qtype) {
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

// --o--s---o--
//   |      |
//   n      n
//   |--G-o-|
//        |
//        s
//        |-discard
// s, n are switches, G is the conditional zx diagram specified by 'left' and
// 'right'. s is on when the input is 0, n is on when the input is 1.
// return the input and the output of the conditional zx along with a vector of
// control spiders
std::pair<std::pair<ZXVert, ZXVert>, std::vector<ZXVert>> add_conditional_zx(
    ZXDiagram& zxd, const ZXVert& left, const ZXVert& right,
    const QuantumType& qtype) {
  ZXVert in = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert out = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert discard_ctr = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  ZXVert discard = zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
  auto switch_s0 = add_switch(zxd, false, qtype);
  auto switch_n0 = add_switch(zxd, true, qtype);
  auto switch_n1 = add_switch(zxd, true, qtype);
  auto switch_s1 = add_switch(zxd, false, qtype);
  zxd.add_wire(in, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out, switch_s0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(in, switch_n0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(left, switch_n0.second, ZXWireType::Basic, qtype);
  zxd.add_wire(right, discard_ctr, ZXWireType::Basic, qtype);
  zxd.add_wire(discard_ctr, switch_n1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(out, switch_n1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(discard_ctr, switch_s1.second, ZXWireType::Basic, qtype);
  zxd.add_wire(discard, switch_s1.second, ZXWireType::Basic, qtype);
  return {
      {in, out},
      {switch_s0.first, switch_s1.first, switch_n0.first, switch_n1.first}};
}

// n-bit AND spider to control the switches
// https://arxiv.org/abs/1910.06818
std::pair<ZXVertPortVec, ZXVertPort> add_n_bit_and(
    ZXDiagram& zxd, unsigned n, const QuantumType& qtype) {
  TKET_ASSERT(n > 1);
  ZXVert z_vert = zxd.add_vertex(ZXType::ZSpider, 0, qtype);
  // Add Triangle -1
  ZXVert z_pi_0 = zxd.add_vertex(ZXType::ZSpider, 1, qtype);
  ZXVert tri_1 = zxd.add_vertex(ZXType::Triangle, qtype);
  ZXVert z_pi_1 = zxd.add_vertex(ZXType::ZSpider, 1, qtype);
  zxd.add_wire(z_vert, z_pi_0, ZXWireType::Basic, qtype);
  zxd.add_wire(tri_1, z_pi_0, ZXWireType::Basic, qtype, 0);
  zxd.add_wire(tri_1, z_pi_1, ZXWireType::Basic, qtype, 1);
  ZXVertPortVec inputs;
  for (unsigned i = 0; i < n; i++) {
    ZXVert tri_0 = zxd.add_vertex(ZXType::Triangle, qtype);
    zxd.add_wire(tri_0, z_vert, ZXWireType::Basic, qtype, 1);
    inputs.push_back({tri_0, 0});
  }
  return {inputs, {z_pi_1, std::nullopt}};
}
// Add converted circ into zxd. Set add_boundary to true
// if boundary spiders are to be added to the boundary.
BoundaryVertMap circuit_to_zx_recursive(
    const Circuit& circ, ZXDiagram& zxd, bool add_boundary) {
  std::map<TypedVertPort, ZXVertPort> vert_lookup;
  std::map<VertPort, ZXVertPort> boolean_outport_lookup;
  BoundaryVertMap bmap;

  // Convert each vertex to ZXDiagram, raise error if not supported
  BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
    // We currently throw an error if the vertex is either classical or flow
    Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
    if (is_flowop_type(op->get_type()) || is_classical_type(op->get_type())) {
      throw Unsupported(
          "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
    switch (op->get_type()) {
      case OpType::Input: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Input, QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::Output: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Output, QuantumType::Quantum);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClInput: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Input, QuantumType::Classical);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        bmap.insert({zx_vert, vert});
        break;
      }
      case OpType::ClOutput: {
        ZXVert zx_vert = zxd.add_vertex(ZXType::Output, QuantumType::Classical);
        if (add_boundary) zxd.add_boundary(zx_vert);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
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
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        zxd.multiply_scalar(0.5);
        break;
      }
      case OpType::Rz: {
        ZXVert zx_vert = zxd.add_vertex(
            ZXType::ZSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        break;
      }
      case OpType::Rx: {
        ZXVert zx_vert = zxd.add_vertex(
            ZXType::XSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        break;
      }
      case OpType::X: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Quantum);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        break;
      }
      case OpType::Z: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::ZSpider, 1, QuantumType::Quantum);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        break;
      }
      case OpType::CX: {
        ZXVert zx_x_vert =
            zxd.add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        ZXVert zx_z_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        zxd.add_wire(zx_x_vert, zx_z_vert);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_z_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_z_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::In}, {zx_x_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::Out}, {zx_x_vert, std::nullopt}});
        zxd.multiply_scalar(2);
        break;
      }
      case OpType::CZ: {
        ZXVert zx_za_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        ZXVert zx_zb_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        zxd.add_wire(zx_za_vert, zx_zb_vert, ZXWireType::H);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_za_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_za_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::In}, {zx_zb_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::Out}, {zx_zb_vert, std::nullopt}});
        break;
      }
      case OpType::Measure: {
        // Add a decoherence node
        ZXVert zx_measure_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        // Add a delete operator
        ZXVert zx_delete_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_measure_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_measure_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::In}, {zx_delete_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 1}, PortType::Out}, {zx_measure_vert, std::nullopt}});
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
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_discard_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_reset_vert, std::nullopt}});
        break;
      }
      case OpType::Collapse: {
        ZXVert zx_vert =
            zxd.add_vertex(ZXType::ZSpider, QuantumType::Classical);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_vert, std::nullopt}});
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_vert, std::nullopt}});
        break;
      }
      case OpType::Create: {
        ZXVert zx_init_vert =
            zxd.add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        zxd.multiply_scalar(0.5);
        vert_lookup.insert(
            {{{vert, 0}, PortType::Out}, {zx_init_vert, std::nullopt}});
        break;
      }
      case OpType::Discard: {
        ZXVert zx_discard_vert =
            zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        vert_lookup.insert(
            {{{vert, 0}, PortType::In}, {zx_discard_vert, std::nullopt}});
        break;
      }
      default:
        if (op->get_type() == OpType::Conditional) {
          // Stores the condition for each port index
          std::vector<bool> port_conditions;

          op_signature_t parent_sig = op->get_signature();

          // Find the nested quantum op and its overall conditions
          // Assume bool ports are always the first few ports
          while (op->get_type() == OpType::Conditional) {
            const Conditional& cond = static_cast<const Conditional&>(*op);
            for (unsigned i = 0; i < cond.get_width(); i++) {
              bool set = cond.get_value() & (1 << (cond.get_width() - i - 1));
              port_conditions.push_back(set);
            }
            op = cond.get_op();
          }
          // Convert the underlying op to zx.
          // If the op is a quantum gate, construct a 1 gate circuit
          // If the op is a box, obtain its circuit decomposition
          op_signature_t inner_sig = op->get_signature();
          Circuit replacement;
          if (is_box_type(op->get_type())) {
            const Box& b = static_cast<const Box&>(*op);
            replacement = *b.to_circuit();
          } else {
            unit_vector_t args;
            unsigned q_index = 0;
            unsigned c_index = 0;
            for (const EdgeType& etype : inner_sig) {
              if (etype == EdgeType::Quantum) {
                args.push_back(Qubit(q_index++));
              } else if (etype == EdgeType::Classical) {
                args.push_back(Bit(c_index++));
              }
              TKET_ASSERT(etype != EdgeType::Boolean);
            }
            replacement = Circuit(q_index, c_index);
            replacement.add_op(op, args);
          }
          BoundaryVertMap box_bm =
              circuit_to_zx_recursive(replacement, zxd, false);

          // The Z spider controlling all switches
          ZXVert master_switch =
              zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);

          // For each qubit/bit path in the inner op, convert it into a
          // conditional path.
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
            // switch_ctrls.insert(switch_ctrls.end(), controls.begin(),
            // controls.end());
            for (const ZXVert& ctr : controls) {
              if (zxd.get_zxtype(ctr) == ZXType::Triangle) {
                zxd.add_wire(ctr, master_switch, ZXWireType::Basic, qtype, 0);
              } else {
                zxd.add_wire(
                    ctr, master_switch, ZXWireType::Basic,
                    QuantumType::Classical);
              }
            }
            // Update lookup
            TKET_ASSERT(parent_sig[i + port_conditions.size()] == inner_sig[i]);
            // Since we assume that the conditions always use the first few
            // ports, the gate path should use the i + port_conditions.size()
            // port.
            vert_lookup.insert(
                {{{vert, i + port_conditions.size()}, PortType::In},
                 {boundary.first, std::nullopt}});
            vert_lookup.insert(
                {{{vert, i + port_conditions.size()}, PortType::Out},
                 {boundary.second, std::nullopt}});
          }
          // Use either a Z spider or a AND operator to connect the master
          // switch and the boolean inputs
          ZXVertPortVec and_inputs;
          if (port_conditions.size() == 1) {
            and_inputs.push_back({master_switch, std::nullopt});
          } else {
            ZXVertPort and_output;
            std::tie(and_inputs, and_output) = add_n_bit_and(
                zxd, port_conditions.size(), QuantumType::Classical);
            zxd.add_wire(
                and_output.first, master_switch, ZXWireType::Basic,
                QuantumType::Classical, and_output.second);
          }
          // Connect inputs to the nodes obtained from above
          for (unsigned i = 0; i < port_conditions.size(); i++) {
            // Each boolean edge shares a source port with other
            // classical/boolean edges. Use the classical Z spider to explicitly
            // implement this copy operation
            ZXVert copy_vert =
                zxd.add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
            Edge in_edge = circ.get_nth_in_edge(vert, i);
            port_t source_port = circ.get_source_port(in_edge);
            Vertex vert_s = circ.source(in_edge);
            // During the wiring stage, each edge originated from {vert_s,
            // p_s} needs go through copy_vert
            boolean_outport_lookup.insert(
                {{vert_s, source_port}, {copy_vert, std::nullopt}});

            if (port_conditions[i]) {
              vert_lookup.insert({{{vert, i}, PortType::In}, and_inputs[i]});
            } else {
              ZXVert x_vert =
                  zxd.add_vertex(ZXType::XSpider, 1, QuantumType::Classical);
              zxd.add_wire(
                  and_inputs[i].first, x_vert, ZXWireType::Basic,
                  QuantumType::Classical, and_inputs[i].second);
              vert_lookup.insert(
                  {{{vert, i}, PortType::In}, {x_vert, std::nullopt}});
            }
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
          // TODO:: We should probably assume that a box can also have boolean
          // input wires op->get_type() == OpType::ClassicalExpBox
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
                {{{vert, port}, PortType::In},
                 {box_bm.right.find(inp)->second, std::nullopt}});
          }
          for (unsigned i = 0; i < q_out_holes.size(); i++) {
            port_t port = circ.get_source_port(q_out_holes[i]);
            Vertex outp = replacement.get_out(Qubit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::Out},
                 {box_bm.right.find(outp)->second, std::nullopt}});
          }
          for (unsigned i = 0; i < c_in_holes.size(); i++) {
            port_t port = circ.get_target_port(c_in_holes[i]);
            Vertex inp = replacement.get_in(Bit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::In},
                 {box_bm.right.find(inp)->second, std::nullopt}});
          }
          for (unsigned i = 0; i < c_out_holes.size(); i++) {
            port_t port = circ.get_source_port(c_out_holes[i]);
            Vertex outp = replacement.get_out(Bit(i));
            vert_lookup.insert(
                {{{vert, port}, PortType::Out},
                 {box_bm.right.find(outp)->second, std::nullopt}});
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
      zxd.add_wire(
          it_s->second.first, it_t->second.first, ZXWireType::Basic,
          QuantumType::Quantum, it_s->second.second, it_t->second.second);
    } else {
      auto bool_it = boolean_outport_lookup.find(VertPort(v_s, p_s));
      if (bool_it == boolean_outport_lookup.end()) {
        zxd.add_wire(
            it_s->second.first, it_t->second.first, ZXWireType::Basic,
            QuantumType::Classical, it_s->second.second, it_t->second.second);
      } else {
        // If the source port is boolean, then connect the copy spider to the
        // source vertex. All out-edges originated from the source port should
        // be connected to the copy spider.
        if (zxd.degree(bool_it->second.first) == 0) {
          zxd.add_wire(
              it_s->second.first, bool_it->second.first, ZXWireType::Basic,
              QuantumType::Classical, it_s->second.second,
              bool_it->second.second);
        }
        zxd.add_wire(
            bool_it->second.first, it_t->second.first, ZXWireType::Basic,
            QuantumType::Classical, bool_it->second.second,
            it_t->second.second);
      }
    }
  }
  return bmap;
}

std::pair<ZXDiagram, BoundaryVertMap> circuit_to_zx(const Circuit& circ) {
  ZXDiagram zxd;
  BoundaryVertMap bmap = circuit_to_zx_recursive(circ, zxd, true);
  // Remove internal boundary vertices produced by the recursion
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
  return {std::move(zxd), std::move(bmap)};
}

}  // namespace tket
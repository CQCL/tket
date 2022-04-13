#include "Circuit/CircPool.hpp"
#include "ZX/Rewrite.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

// enum class PortType { In, Out };
// typedef std::pair<VertPort, PortType> TypedVertPort;

// bool is_spiderless_optype(const OpType& optype) {
//   return optype == OpType::Barrier || optype == OpType::SWAP ||
//          optype == OpType::noop;
// }

// // Construct a switch
// // reverse: true indicates off if 1
// zx::ZXDiagram switch_box(bool reverse) {
//   zx::ZXDiagram zxd(1,2,0,0);
//   zx::ZXVertVec ins = get_boundary(ZXType::Input);
//   zx::ZXVertVec outs = get_boundary(ZXType::Output);
//   zx::ZXVert triangle =
//             zxd.add_vertex(zx::ZXType::Triangle, zx::QuantumType::Quantum);
//   zx::ZXVert x = zxd.add_vertex(zx::ZXType::XSpider, 0,
//   zx::QuantumType::Quantum); zxd.add_wire(triangle, x); zxd.add_wire(outs[0],
//   x); zxd.add_wire(outs[1], x); if (reverse) {
//     zx::ZXVert flip = zxd.add_vertex(zx::ZXType::XSpider, 1,
//     zx::QuantumType::Quantum); zxd.add_wire(triangle, flip);
//     zxd.add_wire(ins[0], flip);
//   }
//   else {
//     zxd.add_wire(ins[0], triangle);
//   }
//   return zxd;
// }

// zx::ZXDiagram conditional_swap(unsigned nwires) {

// }

typedef std::map<
    unsigned, std::pair<std::optional<unsigned>, std::optional<unsigned>>>
    PortMap;


// Convert an Op to a ZXGen
std::tuple<zx::ZXGen_ptr, Expr, std::optional<PortMap>> vertx_to_zx_gen(
    const Op_ptr& op_ptr) {
  // We currently throw an error if the vertex is either
  // 1. box , conditional, classical, flow
  if (is_box_type(op_ptr->get_type()) || is_flowop_type(op_ptr->get_type()) ||
      is_classical_type(op_ptr->get_type()) ||
      op_ptr->get_type() == OpType::Conditional) {
    throw Unsupported(
        "Cannot convert OpType: " + op_ptr->get_name() + " to a ZX node. \n");
  }
  switch (op_ptr->get_type()) {
    case OpType::Input: {
      return {zx::ZXGen::create_gen(zx::ZXType::Input), 1, std::nullopt};
    }
    case OpType::Output: {
      return {zx::ZXGen::create_gen(zx::ZXType::Output), 1, std::nullopt};
    }
    case OpType::ClInput: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::Input, zx::QuantumType::Classical),
          1, std::nullopt};
    }
    case OpType::ClOutput: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::Output, zx::QuantumType::Classical),
          1, std::nullopt};
    }
    case OpType::Barrier: {
      zx::ZXDiagram box;
      PortMap pm;
      op_signature_t sig = op_ptr->get_signature();
      // Add identity for each port
      for (unsigned i = 0; i < sig.size(); i++) {
        zx::QuantumType zx_etype = zx::QuantumType::Classical;
        if (sig[i] == EdgeType::Quantum) {
          zx_etype = zx::QuantumType::Quantum;
        }
        zx::ZXVert in = box.add_vertex(zx::ZXType::Input, zx_etype);
        box.add_boundary(in);
        zx::ZXVert out = box.add_vertex(zx::ZXType::Output, zx_etype);
        box.add_boundary(out);
        box.add_wire(in, out);
        pm.insert({i, {2 * i, 2 * i + 1}});
      }
      return {std::make_shared<const zx::ZXBox>(box), 1, pm};
    }
    case OpType::noop: {
      zx::ZXDiagram box(1, 1, 0, 0);
      PortMap pm;
      pm.insert({0, {0, 1}});
      return {std::make_shared<const zx::ZXBox>(box), 1, std::nullopt};
    }
    case OpType::SWAP: {
      zx::ZXDiagram box(2, 2, 0, 0);
      PortMap pm;
      zx::ZXVertVec ins = box.get_boundary(zx::ZXType::Input);
      zx::ZXVertVec outs = box.get_boundary(zx::ZXType::Output);
      box.add_wire(ins[0], outs[1]);
      box.add_wire(ins[1], outs[0]);
      pm.insert({0, {0, 2}});
      pm.insert({1, {1, 3}});
      return {std::make_shared<const zx::ZXBox>(box), 1, pm};
    }
    case OpType::H: {
      return {zx::ZXGen::create_gen(zx::ZXType::Hbox), 0.5, std::nullopt};
    }
    case OpType::Rz: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::ZSpider, op_ptr->get_params()[0]),
          1, std::nullopt};
    }
    case OpType::Rx: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::XSpider, op_ptr->get_params()[0]),
          1, std::nullopt};
    }
    case OpType::X: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::XSpider, Expr(1)), 1, std::nullopt};
    }
    case OpType::Z: {
      return {
          zx::ZXGen::create_gen(zx::ZXType::ZSpider, Expr(1)), 1, std::nullopt};
    }
    case OpType::CX: {
      zx::ZXDiagram box(2, 2, 0, 0);
      PortMap pm;
      zx::ZXVert zx_z_vert = box.add_vertex(zx::ZXType::ZSpider, 0);
      zx::ZXVert zx_x_vert = box.add_vertex(zx::ZXType::XSpider, 0);
      box.add_wire(zx_x_vert, zx_z_vert);
      zx::ZXVertVec ins = box.get_boundary(zx::ZXType::Input);
      zx::ZXVertVec outs = box.get_boundary(zx::ZXType::Output);
      box.add_wire(ins[0], zx_z_vert);
      box.add_wire(zx_z_vert, outs[0]);
      box.add_wire(ins[1], zx_x_vert);
      box.add_wire(zx_x_vert, outs[1]);
      pm.insert({0, {0, 2}});
      pm.insert({1, {1, 3}});
      return {std::make_shared<const zx::ZXBox>(box), 2, pm};
    }
    case OpType::CZ: {
      zx::ZXDiagram box(2, 2, 0, 0);
      PortMap pm;
      zx::ZXVert zx_za_vert = box.add_vertex(zx::ZXType::ZSpider, 0);
      zx::ZXVert zx_zb_vert = box.add_vertex(zx::ZXType::ZSpider, 0);
      box.add_wire(zx_za_vert, zx_zb_vert, zx::ZXWireType::H);
      zx::ZXVertVec ins = box.get_boundary(zx::ZXType::Input);
      zx::ZXVertVec outs = box.get_boundary(zx::ZXType::Output);
      box.add_wire(ins[0], zx_za_vert);
      box.add_wire(zx_za_vert, outs[0]);
      box.add_wire(ins[1], zx_zb_vert);
      box.add_wire(zx_zb_vert, outs[1]);
      pm.insert({0, {0, 2}});
      pm.insert({1, {1, 3}});
      return {std::make_shared<const zx::ZXBox>(box), 1, pm};
    }
    case OpType::Measure: {
      zx::ZXDiagram box(1, 1, 1, 1);
      PortMap pm;
      // Add a decoherence node
      zx::ZXVert zx_measure_vert =
          box.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
      // Add a delete operator
      zx::ZXVert zx_delete_vert =
          box.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
      zx::ZXVertVec ins = box.get_boundary(zx::ZXType::Input);
      zx::ZXVertVec outs = box.get_boundary(zx::ZXType::Output);
      // The zx constructor add quantum boundaries first
      box.add_wire(ins[0], zx_measure_vert);
      box.add_wire(zx_measure_vert, outs[0]);
      box.add_wire(zx_measure_vert, outs[1]);
      box.add_wire(ins[1], zx_delete_vert);
      pm.insert({0, {0, 1}});
      pm.insert({1, {2, 3}});
      return {std::make_shared<const zx::ZXBox>(box), 1, pm};
    }
    case OpType::Reset: {
      zx::ZXDiagram box(1, 1, 0, 0);
      PortMap pm;
      // Discard
      zx::ZXVert zx_discard_vert =
          box.add_vertex(zx::ZXType::ZSpider, 0, zx::QuantumType::Classical);
      // Add a node to prepare |0>
      zx::ZXVert zx_reset_vert = box.add_vertex(zx::ZXType::XSpider, 0);
      zx::ZXVertVec ins = box.get_boundary(zx::ZXType::Input);
      zx::ZXVertVec outs = box.get_boundary(zx::ZXType::Output);
      box.add_wire(ins[0], zx_discard_vert);
      box.add_wire(zx_reset_vert, outs[0]);
      pm.insert({0, {0, 1}});
      return {std::make_shared<const zx::ZXBox>(box), 0.5, pm};
    }
    case OpType::Collapse: {
      return {
          zx::ZXGen::create_gen(
              zx::ZXType::ZSpider, zx::QuantumType::Classical),
          1, std::nullopt};
    }
    case OpType::Create: {
      return {zx::ZXGen::create_gen(zx::ZXType::XSpider), 0.5, std::nullopt};
    }
    case OpType::Discard: {
      return {
          zx::ZXGen::create_gen(
              zx::ZXType::ZSpider, zx::QuantumType::Classical),
          1, std::nullopt};
    }
    default:
      throw Unsupported(
          "Cannot convert gate type: " + op_ptr->get_name() +
          " to a ZX node, try rebase the gates to use Rx, Rz, X, Z, H, CZ "
          "or CX. \n");
  }
}

zx::ZXDiagram circuit_to_zx(const Circuit& circ) {
  // TODO: how can users know the boundary mapping?
  zx::ZXDiagram zxd;

  sequenced_map_t<Vertex, std::pair<zx::ZXVert, std::optional<PortMap>>>
      vert_lookup;

  // Convert each vertex to ZXDiagram, raise error if not supported
  BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
    Op_ptr op_ptr = circ.get_Op_ptr_from_Vertex(vert);
    zx::ZXGen_ptr zx_gen;
    Expr scalar;
    std::optional<PortMap> pm;
    // They order of the boundary generators matches the ports of the op
    std::tie(zx_gen, scalar, pm) = vertx_to_zx_gen(op_ptr);
    zx::ZXVert zx_vert = zxd.add_vertex(zx_gen);
    if (is_boundary_q_type(op_ptr->get_type()) ||
        is_boundary_c_type(op_ptr->get_type())) {
      zxd.add_boundary(zx_vert);
    }
    zxd.multiply_scalar(scalar);
    vert_lookup.insert({vert, {zx_vert, pm}});
  }

  // Connect the ZX nodes
  BGL_FORALL_EDGES(edge, circ.dag, DAG) {
    Vertex v_s = circ.source(edge);
    Vertex v_t = circ.target(edge);
    port_t p_s = circ.get_source_port(edge);
    port_t p_t = circ.get_target_port(edge);
    std::optional<unsigned> zx_p_s = std::nullopt;
    std::optional<unsigned> zx_p_t = std::nullopt;

    auto it_s = vert_lookup.get<TagKey>().find(v_s);
    auto it_t = vert_lookup.get<TagKey>().find(v_t);

    if (it_s->second.second) {
      zx_p_s = it_s->second.second.value().find(p_s)->second.second;
    }
    if (it_t->second.second) {
      zx_p_t = it_t->second.second.value().find(p_t)->second.first;
    }
    zx::QuantumType qtype = zx::QuantumType::Classical;
    if (circ.get_edgetype(edge) == EdgeType::Quantum) {
      qtype = zx::QuantumType::Quantum;
    }
    zxd.add_wire(
        it_s->second.first, it_t->second.first, zx::ZXWireType::Basic, qtype,
        zx_p_s, zx_p_t);
  }
  // zx::Rewrite::decompose_boxes().apply(zxd);
  return zxd;
}

}  // namespace tket
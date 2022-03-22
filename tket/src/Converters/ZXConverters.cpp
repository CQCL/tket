#include "Circuit/CircPool.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

enum class PortType { In, Out };
typedef std::pair<VertPort, PortType> TypedVertPort;

zx::ZXDiagram circuit_to_zx(const Circuit& circ) {
  zx::ZXDiagram zxd;
  // Raise the scalar to the power of 2 due to doubling
  zxd.multiply_scalar(exp(2. * i_ * PI * circ.get_phase()));

  sequenced_map_t<TypedVertPort, zx::ZXVert> vert_lookup;
  // Append each vertex to ZXDiagram raise error if not supported

  BGL_FORALL_VERTICES(vert, circ.dag, DAG) {
    // We currently throw an error if the vertex is either
    // 1. box , conditional, classical, flow, projective
    Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
    if (is_box_type(op->get_type()) || is_flowop_type(op->get_type()) ||
        is_classical_type(op->get_type()) ||
        op->get_type() == OpType::Conditional) {
      throw Unsupported(
          "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
    switch (op->get_type()) {
      case OpType::Input: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Input, zx::QuantumType::Quantum);
        zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Output: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Output, zx::QuantumType::Quantum);
        zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        break;
      }
      case OpType::ClInput: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Input, zx::QuantumType::Classical);
        zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::ClOutput: {
        zx::ZXVert zx_vert =
            zxd.add_vertex(zx::ZXType::Output, zx::QuantumType::Classical);
        zxd.add_boundary(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        break;
      }
      case OpType::Barrier:
        continue;
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
        zxd.multiply_scalar(exp(-i_ * PI * op->get_params()[0]));
        break;
      }
      case OpType::Rx: {
        zx::ZXVert zx_vert = zxd.add_vertex(
            zx::ZXType::XSpider, op->get_params()[0], zx::QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        zxd.multiply_scalar(exp(-i_ * PI * op->get_params()[0]));
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
      default:
        throw Unsupported(
            "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
  }
  BGL_FORALL_EDGES(edge, circ.dag, DAG) {
    Vertex v_s = circ.source(edge);
    Vertex v_t = circ.target(edge);
    port_t p_s = circ.get_source_port(edge);
    port_t p_t = circ.get_target_port(edge);
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
  return zxd;
}

}  // namespace tket
#include "Circuit/CircPool.hpp"
#include "Transformations/Decomposition.hpp"
#include "Transformations/Rebase.hpp"
#include "Transformations/Transform.hpp"
#include "ZX/ZXDiagram.hpp"

namespace tket {

namespace zx {

enum class PortType { In, Out };
typedef std::pair<VertPort, PortType> TypedVertPort;

ZXDiagram::ZXDiagram(const Circuit& circ) : ZXDiagram() {
  // Make a copy of the original circuit
  Circuit c(circ);
  // Rebase circuit to Rz, Rx, Hadamard, CZ, and CNOTs
  OpTypeSet gates = {OpType::Rz, OpType::Rx, OpType::H, OpType::CX, OpType::CZ};
  Transform t_decompose_box = Transforms::decomp_boxes();
  Transform t_rebase =
      Transforms::rebase_factory(gates, CircPool::CX(), CircPool::tk1_to_rzrx);
  t_decompose_box.apply(c);
  t_rebase.apply(c);
  multiply_scalar(exp(i_ * PI * c.get_phase()));

  sequenced_map_t<TypedVertPort, ZXVert> vert_lookup;
  // Append each vertex to ZXDiagram raise error if not supported

  BGL_FORALL_VERTICES(vert, c.dag, DAG) {
    // We currently throw an error if the vertex is either
    // 1. box , conditional, classical, flow, projective
    Op_ptr op = c.get_Op_ptr_from_Vertex(vert);
    if (is_box_type(op->get_type()) ||
        is_flowop_type(op->get_type()) || is_classical_type(op->get_type()) ||
        op->get_type() == OpType::Conditional) {
      throw Unsupported(
          "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
    switch (op->get_type()) {
      case OpType::Input: {
        ZXVert zx_vert = add_vertex(ZXType::Input, QuantumType::Quantum);
        boundary.push_back(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Output: {
        ZXVert zx_vert = add_vertex(ZXType::Output, QuantumType::Quantum);
        boundary.push_back(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        break;
      }
      case OpType::ClInput: {
        ZXVert zx_vert = add_vertex(ZXType::Input, QuantumType::Classical);
        boundary.push_back(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::ClOutput: {
        ZXVert zx_vert = add_vertex(ZXType::Output, QuantumType::Classical);
        boundary.push_back(zx_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        break;
      }
      case OpType::Barrier:
        continue;
      case OpType::H: {
        ZXVert zx_vert = add_vertex(ZXType::Hbox, QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      case OpType::Rz: {
        ZXVert zx_vert = add_vertex(
            ZXType::ZSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        multiply_scalar(exp(-i_ * 0.5* PI * op->get_params()[0]));
        break;
      }
      case OpType::Rx: {
        ZXVert zx_vert = add_vertex(
            ZXType::XSpider, op->get_params()[0], QuantumType::Quantum);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        multiply_scalar(exp(-i_ * 0.5* PI * op->get_params()[0]));
        break;
      }
      case OpType::CX: {
        ZXVert zx_x_vert = add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        ZXVert zx_z_vert = add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        add_wire(zx_x_vert, zx_z_vert);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_z_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_z_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_x_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_x_vert});
        multiply_scalar(sqrt(2));
        break;
      }
      case OpType::CZ: {
        ZXVert zx_za_vert =
            add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        ZXVert zx_zb_vert =
            add_vertex(ZXType::ZSpider, 0, QuantumType::Quantum);
        add_wire(zx_za_vert, zx_zb_vert, ZXWireType::H);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_za_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_za_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_zb_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_zb_vert});
        multiply_scalar(sqrt(2));
        break;
      }
      case OpType::Measure: {
        // Add a decoherence node
        ZXVert zx_measure_vert =
            add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        // Add a delete operator
        ZXVert zx_delete_vert = add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_measure_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_measure_vert});
        vert_lookup.insert({{{vert, 1}, PortType::In}, zx_delete_vert});
        vert_lookup.insert({{{vert, 1}, PortType::Out}, zx_measure_vert});
        break;
      }
      case OpType::Reset: {
        // Discard
        ZXVert zx_discard_vert =
            add_vertex(ZXType::ZSpider, 0, QuantumType::Classical);
        // Add a node to prepare |0>
        ZXVert zx_reset_vert = add_vertex(ZXType::XSpider, 0, QuantumType::Quantum);
        multiply_scalar(sqrt(0.5));
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_discard_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_reset_vert});
        break;
      }
      case OpType::Collapse: {
        ZXVert zx_vert = add_vertex(ZXType::ZSpider, QuantumType::Classical);
        vert_lookup.insert({{{vert, 0}, PortType::In}, zx_vert});
        vert_lookup.insert({{{vert, 0}, PortType::Out}, zx_vert});
        break;
      }
      default:
        throw Unsupported(
            "Cannot convert OpType: " + op->get_name() + " to a ZX node. \n");
    }
  }
  BGL_FORALL_EDGES(edge, c.dag, DAG) {
    Vertex v_s = c.source(edge);
    Vertex v_t = c.target(edge);
    port_t p_s = c.get_source_port(edge);
    port_t p_t = c.get_target_port(edge);
    auto it_s = vert_lookup.get<TagKey>().find(TypedVertPort(VertPort(v_s, p_s), PortType::Out));
    auto it_t = vert_lookup.get<TagKey>().find(TypedVertPort(VertPort(v_t, p_t), PortType::In));
    if(c.get_edgetype(edge) == EdgeType::Quantum) {
      add_wire(it_s->second, it_t->second);
    }
    else {
      add_wire(it_s->second, it_t->second, ZXWireType::Basic, QuantumType::Classical);
    }
    
  }
}

}  // namespace zx

}  // namespace tket
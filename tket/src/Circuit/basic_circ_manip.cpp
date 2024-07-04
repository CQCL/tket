// Copyright 2019-2024 Cambridge Quantum Computing
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

/////////////////////////////////////////////////////
// ALL METHODS TO PERFORM BASIC CIRCUIT MANIPULATION//
/////////////////////////////////////////////////////

#include <memory>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "tket/Circuit/Boxes.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Ops/BarrierOp.hpp"
#include "tket/Ops/MetaOp.hpp"

namespace tket {

// if there are any blank wires in the circuit,
// this method removes them and removes the vertices
// from boundaries
void Circuit::remove_blank_wires(bool keep_blank_classical_wires) {
  bool found_empty_bit_at_end = true;
  while (found_empty_bit_at_end) {
    found_empty_bit_at_end = false;
    VertexList bin;
    unit_vector_t unused_units;
    const Op_ptr noop = get_op_ptr(OpType::noop);

    for (const BoundaryElement& el : boundary.get<TagID>()) {
      if (!keep_blank_classical_wires || el.type() == UnitType::Qubit) {
        Vertex in = el.in_;
        Vertex out = el.out_;
        VertexVec succs = get_successors(in);

        bool remove_op = false;

        // check if the unitid is unused
        if (succs.front() == out && succs.size() == 1) {
          remove_op = true;
        }

        // check if the unused unitid is a bit in a register of dim 1
        if (remove_op && el.type() == UnitType::Bit && el.id_.reg_dim() == 1) {
          // check if the empty bit is at the end of a register
          if (get_reg(el.id_.reg_name()).size() != (el.id_.index()[0] + 1)) {
            remove_op = false;
          } else {
            // rerun while loop
            found_empty_bit_at_end = true;
          }
        }

        if (remove_op) {
          dag[in].op = noop;
          bin.push_back(in);
          dag[out].op = noop;
          bin.push_back(out);
          unused_units.push_back(el.id_);
        }
      }
    }
    for (const UnitID& u : unused_units) {
      boundary.get<TagID>().erase(u);
    }
    remove_vertices(bin, GraphRewiring::No, VertexDeletion::Yes);
  }
}

void Circuit::remove_noops() {
  VertexSet bin;
  BGL_FORALL_VERTICES(v, dag, DAG) {
    Op_ptr op = get_Op_ptr_from_Vertex(v);
    std::optional<double> phase;
    if (op->get_desc().is_gate() && (phase = op->is_identity())) {
      remove_vertex(
          v, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
      add_phase(*phase);
      bin.insert(v);
    }
  }
  remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
}

std::ostream& operator<<(std::ostream& out, const Circuit& circ) {
  for (const Command& com : circ) out << com << std::endl;
  out << "Phase (in half-turns): " << circ.get_phase() << std::endl;
  return out;
}

template <>
Vertex Circuit::add_op<unsigned>(
    const Op_ptr& gate, const std::vector<unsigned>& args,
    std::optional<std::string> opgroup) {
  op_signature_t sig = gate->get_signature();

  if (sig.size() != args.size()) {
    throw CircuitInvalidity(
        std::to_string(args.size()) + " args provided, but " +
        gate->get_name() + " requires " + std::to_string(sig.size()));
  }

  OpType optype = gate->get_type();
  unit_vector_t arg_ids;
  for (unsigned i = 0; i < args.size(); ++i) {
    switch (sig.at(i)) {
      case EdgeType::Quantum: {
        arg_ids.push_back(Qubit(args[i]));
        break;
      }
      case EdgeType::Classical:
      case EdgeType::Boolean: {
        arg_ids.push_back(Bit(args[i]));
        break;
      }
      default: {
        TKET_ASSERT(!"add_op found invalid edge type in signature");
      }
    }
  }
  if (optype == OpType::CnRy && args.size() == 1) {
    return add_op(get_op_ptr(OpType::Ry, gate->get_params()), arg_ids);
  } else if (optype == OpType::CnRx && args.size() == 1) {
    return add_op(get_op_ptr(OpType::Rx, gate->get_params()), arg_ids);
  } else if (optype == OpType::CnRz && args.size() == 1) {
    return add_op(get_op_ptr(OpType::Rz, gate->get_params()), arg_ids);
  } else if (optype == OpType::CnX && args.size() == 1) {
    return add_op(get_op_ptr(OpType::X), arg_ids);
  } else if (optype == OpType::CnZ && args.size() == 1) {
    return add_op(get_op_ptr(OpType::Z), arg_ids);
  } else if (optype == OpType::CnY && args.size() == 1) {
    return add_op(get_op_ptr(OpType::Y), arg_ids);
  }
  return add_op(gate, arg_ids, opgroup);
}

Vertex Circuit::add_barrier(
    const std::vector<unsigned>& qubits, const std::vector<unsigned>& bits,
    const std::string& _data) {
  op_signature_t sig(qubits.size(), EdgeType::Quantum);
  op_signature_t cl_sig(bits.size(), EdgeType::Classical);
  sig.insert(sig.end(), cl_sig.begin(), cl_sig.end());
  std::vector<unsigned> args = qubits;
  args.insert(args.end(), bits.begin(), bits.end());
  return add_op(std::make_shared<BarrierOp>(sig, _data), args);
}

Vertex Circuit::add_barrier(
    const unit_vector_t& args, const std::string& _data) {
  op_signature_t sig;
  for (const UnitID& arg : args) {
    if (arg.type() == UnitType::Qubit) {
      sig.push_back(EdgeType::Quantum);
    } else {
      sig.push_back(EdgeType::Classical);
    }
  }
  return add_op(std::make_shared<BarrierOp>(sig, _data), args);
}

Vertex Circuit::add_conditional_barrier(
    const std::vector<unsigned>& barrier_qubits,
    const std::vector<unsigned>& barrier_bits,
    const std::vector<unsigned>& condition_bits, unsigned value,
    const std::string& _data, std::optional<std::string> opgroup) {
  op_signature_t sig(barrier_qubits.size(), EdgeType::Quantum);
  sig.insert(sig.end(), barrier_bits.size(), EdgeType::Classical);

  std::vector<unsigned> args = condition_bits;
  args.insert(args.end(), barrier_qubits.begin(), barrier_qubits.end());
  args.insert(args.end(), barrier_bits.begin(), barrier_bits.end());
  return add_op(
      std::make_shared<Conditional>(
          std::make_shared<BarrierOp>(sig, _data),
          (unsigned)condition_bits.size(), value),
      args, opgroup);
}

Vertex Circuit::add_conditional_barrier(
    const unit_vector_t& barrier_args, const bit_vector_t& condition_bits,
    unsigned value, const std::string& _data,
    std::optional<std::string> opgroup) {
  op_signature_t sig;
  for (const UnitID& arg : barrier_args) {
    if (arg.type() == UnitType::Qubit) {
      sig.push_back(EdgeType::Quantum);
    } else {
      TKET_ASSERT(arg.type() == UnitType::Bit);
      sig.push_back(EdgeType::Classical);
    }
  }
  unit_vector_t args;
  args.insert(args.end(), condition_bits.begin(), condition_bits.end());
  args.insert(args.end(), barrier_args.begin(), barrier_args.end());
  return add_op(
      std::make_shared<Conditional>(
          std::make_shared<BarrierOp>(sig, _data),
          (unsigned)condition_bits.size(), value),
      args, opgroup);
}

std::string Circuit::get_next_c_reg_name(const std::string& reg_name) {
  if (!get_reg_info(reg_name)) {
    return reg_name;
  }
  unsigned post_fix = 1;
  while (true) {
    std::string incremented_reg_name =
        reg_name + "(" + std::to_string(post_fix) + ")";
    if (!get_reg_info(incremented_reg_name)) {
      return incremented_reg_name;
    }
    post_fix++;
  }
}

static void append_debug_bits(
    Circuit& circ, unit_vector_t& args,
    const std::vector<bool>& expected_readouts,
    const std::optional<std::string>& postfix) {
  unsigned n_temp_bits = expected_readouts.size();
  unsigned n_one_bits =
      std::accumulate(expected_readouts.begin(), expected_readouts.end(), 0);
  unsigned n_zero_bits = n_temp_bits - n_one_bits;

  // Add a classical registers for this assertion
  std::string zero_reg_name =
      (postfix == std::nullopt)
          ? c_debug_zero_prefix() + "_" + c_debug_default_name()
          : c_debug_zero_prefix() + "_" + postfix.value();
  std::string one_reg_name =
      (postfix == std::nullopt)
          ? c_debug_one_prefix() + "_" + c_debug_default_name()
          : c_debug_one_prefix() + "_" + postfix.value();
  std::string incremented_zero_reg_name =
      circ.get_next_c_reg_name(zero_reg_name);
  if (n_zero_bits > 0) {
    circ.add_c_register(incremented_zero_reg_name, n_zero_bits);
  }
  std::string incremented_one_reg_name = circ.get_next_c_reg_name(one_reg_name);
  if (n_one_bits > 0) {
    circ.add_c_register(incremented_one_reg_name, n_one_bits);
  }

  unsigned zero_reg_index = 0;
  unsigned one_reg_index = 0;

  for (unsigned i = 0; i < n_temp_bits; i++) {
    if (expected_readouts[i]) {
      args.push_back(Bit(incremented_one_reg_name, one_reg_index++));
    } else {
      args.push_back(Bit(incremented_zero_reg_name, zero_reg_index++));
    }
  }
}

Vertex Circuit::add_assertion(
    const ProjectorAssertionBox& assertion_box,
    const std::vector<Qubit>& qubits, const std::optional<Qubit>& ancilla,
    const std::optional<std::string>& name) {
  auto circ_ptr = assertion_box.to_circuit();
  unsigned log2_dim = log2(assertion_box.get_matrix().rows());

  if (circ_ptr->n_qubits() > log2_dim && ancilla == std::nullopt) {
    throw CircuitInvalidity("This assertion requires an ancilla");
  }

  if (qubits.size() != log2_dim) {
    throw CircuitInvalidity(
        std::to_string(qubits.size()) +
        " target qubits provided, but the projector requires " +
        std::to_string(log2_dim));
  }

  unit_vector_t args;
  args.insert(args.end(), qubits.begin(), qubits.end());
  if (circ_ptr->n_qubits() > log2_dim) {
    args.push_back(*ancilla);
  }
  append_debug_bits(*this, args, assertion_box.get_expected_readouts(), name);
  return add_op<UnitID>(
      std::make_shared<ProjectorAssertionBox>(assertion_box), args);
}

Vertex Circuit::add_assertion(
    const StabiliserAssertionBox& assertion_box,
    const std::vector<Qubit>& qubits, const Qubit& ancilla,
    const std::optional<std::string>& name) {
  auto circ_ptr = assertion_box.to_circuit();
  unsigned pauli_len = assertion_box.get_stabilisers()[0].string.size();
  if (qubits.size() != pauli_len) {
    throw CircuitInvalidity(
        std::to_string(qubits.size()) +
        " target qubits provided, but the stabilisers requires " +
        std::to_string(pauli_len));
  }

  unit_vector_t args;
  args.insert(args.end(), qubits.begin(), qubits.end());
  args.push_back(ancilla);
  append_debug_bits(*this, args, assertion_box.get_expected_readouts(), name);
  return add_op<UnitID>(
      std::make_shared<StabiliserAssertionBox>(assertion_box), args);
}

// adds a vertex to dag of given op type without any connecting edges
// does not add boundary vertices to registers; this should be done manually
Vertex Circuit::add_vertex(
    const Op_ptr op_ptr, std::optional<std::string> opgroup) {
  Vertex new_V = boost::add_vertex(this->dag);
  this->dag[new_V] = {op_ptr, opgroup};
  return new_V;
}

Vertex Circuit::add_vertex(
    const OpType& type, std::optional<std::string> opgroup) {
  return add_vertex(get_op_ptr(type), opgroup);
}

// given vertices and desired in ports for i2 and out ports for i1, adds
// edge between them
// there are no checks to ensure the vertex exists in the graph
Edge Circuit::add_edge(
    const VertPort& source, const VertPort& target, const EdgeType& type) {
  // add new edge
  std::pair<Edge, bool> edge_pairy =
      boost::add_edge(source.first, target.first, this->dag);

  TKET_ASSERT(edge_pairy.second);  // Cannot create edge between vertices

  Edge new_E = edge_pairy.first;
  dag[new_E].ports.first = source.second;
  dag[new_E].ports.second = target.second;
  dag[new_E].type = type;

  return new_E;
}

// given vertex, eradicates it from dag, returning edges used for rewiring
// methods
// there are no checks to ensure the vertex exists in the graph. This may cause
// a segfault if called on a vertex which doesn't exist any more. Take care with
// ordering of vertex removal! has flags for whether you want to a) rewire the
// hole left by the vertex and b) destroy vertex once it has been isolated As a
// sanity check, you cannot remove boundary vertices; converting them to noops
// beforehand is recommended
void Circuit::remove_vertex(
    const Vertex& deadvert, GraphRewiring graph_rewiring,
    VertexDeletion vertex_deletion) {
  if (graph_rewiring == GraphRewiring::Yes) {
    EdgeVec ins = get_in_edges(deadvert);
    std::vector<EdgeVec> bundles = get_b_out_bundles(deadvert);
    port_t p = 0;
    for (const Edge& e : ins) {
      EdgeType type = get_edgetype(e);
      if (type != EdgeType::Boolean) {
        Vertex pred_vert = source(e);
        port_t pred_port = get_source_port(e);
        Edge next_e = get_nth_out_edge(deadvert, p);
        Vertex succ_vert = target(next_e);
        port_t succ_port = get_target_port(next_e);
        add_edge({pred_vert, pred_port}, {succ_vert, succ_port}, type);
        if (type == EdgeType::Classical) {
          for (const Edge& cr : bundles[p]) {
            Vertex sv = target(cr);
            port_t sp = get_target_port(cr);
            add_edge({pred_vert, pred_port}, {sv, sp}, EdgeType::Boolean);
          }
        }
      }
      ++p;
    }
  }

  boost::clear_vertex(deadvert, this->dag);
  if (vertex_deletion == VertexDeletion::Yes) {
    // Cannot remove a boundary vertex
    TKET_ASSERT(!detect_boundary_Op(deadvert));
    boost::remove_vertex(deadvert, this->dag);
  }
}

// same as previous but for a set of vertices
void Circuit::remove_vertices(
    const VertexSet& surplus, GraphRewiring graph_rewiring,
    VertexDeletion vertex_deletion) {
  for (const Vertex& to_remove : surplus) {
    remove_vertex(to_remove, graph_rewiring, vertex_deletion);
  }
}

void Circuit::remove_vertices(
    const VertexList& surplus, GraphRewiring graph_rewiring,
    VertexDeletion vertex_deletion) {
  for (const Vertex& to_remove : surplus) {
    remove_vertex(to_remove, graph_rewiring, vertex_deletion);
  }
}

void Circuit::remove_edge(const Edge& edge) {
  boost::remove_edge(edge, this->dag);
}

unit_map_t Circuit::flatten_registers() {
  unit_map_t rename_map;
  unsigned q_index = 0;
  unsigned c_index = 0;
  for (const BoundaryElement& el : boundary.get<TagID>()) {
    if (el.type() == UnitType::Qubit) {
      rename_map.insert({el.id_, Qubit(q_index++)});
    } else {
      rename_map.insert({el.id_, Bit(c_index++)});
    }
  }
  try {
    rename_units(rename_map);
  } catch (const std::exception& e) {
    std::stringstream ss;
    ss << "Unable to flatten registers: " << e.what();
    throw std::runtime_error(ss.str());
  }
  return rename_map;
}

// this automatically updates the circuit boundaries
void Circuit::add_blank_wires(unsigned n) {
  TKET_ASSERT(default_regs_ok());  // Incompatible registers exist with the
                                   // default names

  unsigned index = 0;
  for (unsigned i = 0; i < n; i++) {
    Vertex in = add_vertex(OpType::Input);
    Vertex out = add_vertex(OpType::Output);
    add_edge({in, 0}, {out, 0}, EdgeType::Quantum);
    bool searching = true;
    while (searching) {
      Qubit q_id(index);
      boundary_t::iterator found = boundary.get<TagID>().find(q_id);
      if (found == boundary.get<TagID>().end()) {
        searching = false;
        boundary.insert({q_id, in, out});
      }
      index++;
    }
  }
}

void Circuit::add_qubit(const Qubit& id, bool reject_dups) {
  boundary_t::index<TagID>::type::iterator found =
      boundary.get<TagID>().find(id);
  if (found != boundary.get<TagID>().end()) {
    if (reject_dups) {
      throw CircuitInvalidity(
          "A unit with ID \"" + id.repr() + "\" already exists");
    } else if (found->type() == UnitType::Qubit) {
      return;
    }
  }
  opt_reg_info_t reg_info = get_reg_info(id.reg_name());
  register_info_t correct_info = {UnitType::Qubit, id.reg_dim()};

  if (reg_info && !(reg_info.value() == correct_info)) {
    throw CircuitInvalidity(
        "Cannot add qubit with ID \"" + id.repr() +
        "\" as register is not compatible");
  }

  Vertex in = add_vertex(OpType::Input);
  Vertex out = add_vertex(OpType::Output);
  add_edge({in, 0}, {out, 0}, EdgeType::Quantum);
  boundary.insert({id, in, out});
}

void Circuit::add_bit(const Bit& id, bool reject_dups) {
  boundary_t::index<TagID>::type::iterator found =
      boundary.get<TagID>().find(id);
  if (found != boundary.get<TagID>().end()) {
    if (reject_dups) {
      throw CircuitInvalidity(
          "A unit with ID \"" + id.repr() + "\" already exists");
    } else if (found->type() == UnitType::Bit) {
      return;
    }
  }
  opt_reg_info_t reg_info = get_reg_info(id.reg_name());
  register_info_t correct_info = {UnitType::Bit, id.reg_dim()};

  if (reg_info && !(reg_info.value() == correct_info)) {
    throw CircuitInvalidity(
        "Cannot add bit with ID \"" + id.repr() +
        "\" as register is not compatible");
  }

  Vertex in = add_vertex(OpType::ClInput);
  Vertex out = add_vertex(OpType::ClOutput);
  add_edge({in, 0}, {out, 0}, EdgeType::Classical);
  boundary.insert({id, in, out});
}

register_t Circuit::add_q_register(std::string reg_name, unsigned size) {
  if (get_reg_info(reg_name))
    throw CircuitInvalidity(
        "A q register with name \"" + reg_name + "\" already exists");
  register_t ids;
  for (unsigned i = 0; i < size; i++) {
    Vertex in = add_vertex(OpType::Input);
    Vertex out = add_vertex(OpType::Output);
    add_edge({in, 0}, {out, 0}, EdgeType::Quantum);
    Qubit id(reg_name, i);
    boundary.insert({id, in, out});
    ids.insert({i, id});
  }
  return ids;
}

register_t Circuit::add_c_register(std::string reg_name, unsigned size) {
  if (get_reg_info(reg_name))
    throw CircuitInvalidity(
        "A c register with name \"" + reg_name + "\" already exists");
  register_t ids;
  for (unsigned i = 0; i < size; i++) {
    Vertex in = add_vertex(OpType::ClInput);
    Vertex out = add_vertex(OpType::ClOutput);
    add_edge({in, 0}, {out, 0}, EdgeType::Classical);
    Bit id(reg_name, i);
    boundary.insert({id, in, out});
    ids.insert({i, id});
  }
  return ids;
}

void Circuit::add_wasm_register(std::size_t number_of_w_) {
  while (number_of_w_ > _number_of_wasm_wires) {
    Vertex in = add_vertex(OpType::WASMInput);
    Vertex out = add_vertex(OpType::WASMOutput);
    add_edge({in, 0}, {out, 0}, EdgeType::WASM);
    WasmState wuid = WasmState(_number_of_wasm_wires);
    wasmwire.push_back(wuid);
    boundary.insert({wuid, in, out});
    ++_number_of_wasm_wires;
  }
}

void Circuit::qubit_create(const Qubit& id) {
  Vertex v = get_in(id);
  dag[v].op = std::make_shared<const MetaOp>(OpType::Create);
}

void Circuit::qubit_create_all() {
  for (const auto& qb : all_qubits()) {
    qubit_create(qb);
  }
}

void Circuit::qubit_discard(const Qubit& id) {
  Vertex v = get_out(id);
  dag[v].op = std::make_shared<const MetaOp>(OpType::Discard);
}

void Circuit::qubit_discard_all() {
  for (const auto& qb : all_qubits()) {
    qubit_discard(qb);
  }
}

// for wiring in a single vertex with multiple qubits
// there are no checks to ensure the vertex exists in the graph
void Circuit::rewire(
    const Vertex& new_vert, const EdgeVec& preds, const op_signature_t& types) {
  // multi qubit gate
  EdgeList bin;
  for (unsigned i = 0; i < preds.size(); ++i) {
    EdgeType insert_type = types[i];
    EdgeType replace_type = get_edgetype(preds[i]);
    port_t port1 = get_source_port(preds[i]);
    port_t port2 = get_target_port(preds[i]);
    Vertex old_v1 = source(preds[i]);
    Vertex old_v2 = target(preds[i]);

    if (insert_type == EdgeType::Boolean) {
      // Cannot rewire; Boolean needs a classical value to read from
      TKET_ASSERT(replace_type == EdgeType::Classical);

      add_edge({old_v1, port1}, {new_vert, i}, insert_type);
    } else {
      // Cannot rewire; type of edge
      if (insert_type != replace_type) {
        throw CircuitInvalidity(
            "Operation can not be added, found invalid parameter type.");
      }

      add_edge({old_v1, port1}, {new_vert, i}, insert_type);
      add_edge({new_vert, i}, {old_v2, port2}, insert_type);
      bin.push_back(preds[i]);
    }
  }
  for (const Edge& e : bin) remove_edge(e);
}

}  // namespace tket

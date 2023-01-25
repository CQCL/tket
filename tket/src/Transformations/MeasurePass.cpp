// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "MeasurePass.hpp"

#include <optional>
#include <tuple>

#include "Circuit/DAGDefs.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Ops/Op.hpp"
#include "Transform.hpp"

namespace tket {

namespace Transforms {

Transform delay_measures() {
  return Transform(
      [](Circuit& circ) { return DelayMeasures::run_delay_measures(circ, false).first; });
}

namespace DelayMeasures {

bool check_only_end_measures(const Command& com, unit_set_t& measured_units) {
  OpType optype = com.get_op_ptr()->get_type();
  if (optype == OpType::Conditional) {
    unit_vector_t all_args = com.get_args();
    const Conditional& cond =
        static_cast<const Conditional&>(*com.get_op_ptr());
    unit_vector_t::iterator arg_it = all_args.begin();
    for (unsigned i = 0; i < cond.get_width(); ++i) {
      if (measured_units.find(*arg_it) != measured_units.end()) return false;
      ++arg_it;
    }
    unit_vector_t new_args = {arg_it, all_args.end()};
    Command new_com = {cond.get_op(), new_args};
    return check_only_end_measures(new_com, measured_units);
  } else if (optype == OpType::CircBox || optype == OpType::CustomGate) {
    const Box& box = static_cast<const Box&>(*com.get_op_ptr());
    unit_map_t interface;
    unit_set_t inner_set;
    unsigned q_count = 0;
    unsigned b_count = 0;
    for (const UnitID& u : com.get_args()) {
      UnitID inner_unit = (u.type() == UnitType::Qubit)
                              ? static_cast<UnitID>(Qubit(q_count++))
                              : static_cast<UnitID>(Bit(b_count++));
      interface.insert({inner_unit, u});
      if (measured_units.find(u) != measured_units.end()) {
        inner_set.insert(inner_unit);
      }
    }
    for (const Command& c : *box.to_circuit()) {
      if (!check_only_end_measures(c, inner_set)) return false;
    }
    for (const UnitID& u : inner_set) {
      measured_units.insert(interface.at(u));
    }
    return true;
  } else if (optype == OpType::Measure) {
    std::pair<unit_set_t::iterator, bool> q_inserted =
        measured_units.insert(com.get_args().at(0));
    std::pair<unit_set_t::iterator, bool> c_inserted =
        measured_units.insert(com.get_args().at(1));
    return q_inserted.second && c_inserted.second;
  } else {
    for (const UnitID& a : com.get_args()) {
      if (measured_units.find(a) != measured_units.end()) return false;
    }
    return true;
  }
}

static Command command_from_vertex(const Circuit& circ, Vertex v) {
  Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
  unit_vector_t args = circ.vertex_unit_map(v);
  Command new_com = {op, args};
  return new_com;
}

/**
 * Follows and edge in the circuit's DAG until it reaches a non-commuting
 * operation.
 * @param circ The circuit
 * @param current_edge The edge to start following.
 * @param commute_pauliZ Whether to allow gates that commute on the Pauli Z
 *basis.
 * @return The edge leading to the first non-commuting operation.
 **/
static Edge follow_until_noncommuting(
    const Circuit& circ, Edge current_edge, bool commute_pauliZ);

/**
 * Checks if the operation commutes with the Pauli Z basis.
 * @param circ The circuit containing the operation.
 * @param v The vertex of the operation.
 * @param in_port The input port of the operation to check.
 * @param commute_pauliZ Whether to allow gates that commute on the Pauli Z
 *basis.
 * @return The output port of the operation if it commutes, or nullopt.
 **/
static std::optional<port_t> op_commutes(
    const Circuit& circ, const Command& com, port_t in_port,
    bool commute_pauliZ) {
  Op_ptr op = com.get_op_ptr();
  OpType optype = op->get_type();
  UnitID in_unit = com.get_args().at(in_port);
  if (optype == OpType::SWAP) {
    return 1 - in_port;
  } else if (optype == OpType::Conditional) {
    unit_vector_t all_args = com.get_args();
    const Conditional& cond = static_cast<const Conditional&>(*op);
    // Construct the internal command and check it
    // If the port is in the condition, return false
    unit_vector_t::iterator arg_it = all_args.begin();
    for (unsigned i = 0; i < cond.get_width(); ++i) {
      if (*arg_it == in_unit) return std::nullopt;
      ++arg_it;
    }
    unit_vector_t new_args = {arg_it, all_args.end()};
    Command new_com = {cond.get_op(), new_args};
    return op_commutes(
        circ, new_com, in_port - cond.get_width(), commute_pauliZ);
  } else if (optype == OpType::CircBox || optype == OpType::CustomGate) {
    const Box& box = static_cast<const Box&>(*com.get_op_ptr());
    const Circuit& circ = *box.to_circuit();
    Vertex inner_vertex = circ.all_inputs().at(in_port);
    Edge inner_edge = circ.get_nth_out_edge(inner_vertex, 0);
    inner_edge = follow_until_noncommuting(circ, inner_edge, commute_pauliZ);
    inner_vertex = circ.target(inner_edge);
    if (!is_final_q_type(circ.get_OpType_from_Vertex(inner_vertex))) {
      // The operation does not commute through the box
      return std::nullopt;
    }
    auto outputs = circ.all_outputs();
    port_t box_out_port = std::distance(
        outputs.begin(),
        std::find(outputs.begin(), outputs.end(), inner_vertex));
    return box_out_port;
  } else {
    Vertex v = com.get_vertex();
    Pauli basis = commute_pauliZ ? Pauli::Z : Pauli::I;
    if (circ.commutes_with_basis(v, basis, PortType::Target, in_port)) {
      return in_port;
    }
    return std::nullopt;
  }
}

static Edge follow_until_noncommuting(
    const Circuit& circ, Edge current_edge, bool commute_pauliZ) {
  port_t current_port = circ.get_target_port(current_edge);
  Vertex current_vertex = circ.target(current_edge);
  OpType current_optype = circ.get_OpType_from_Vertex(current_vertex);
  while (!is_final_q_type(current_optype)) {
    Command com = command_from_vertex(circ, current_vertex);
    std::optional<port_t> next_port =
        op_commutes(circ, com, current_port, commute_pauliZ);
    if (next_port.has_value()) {
      current_port = next_port.value();
    } else {
      break;
    }
    // Update to successor
    current_edge = circ.get_nth_out_edge(current_vertex, current_port);
    current_vertex = circ.target(current_edge);
    current_port = circ.get_target_port(current_edge);
    current_optype = circ.get_OpType_from_Vertex(current_vertex);
  }
  return current_edge;
}

std::pair<bool, bool> run_delay_measures(Circuit& circ, bool dry_run) {
  bool modified = false;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    OpType optype = circ.get_OpType_from_Vertex(v);
    if (optype == OpType::Measure) {
      Edge c_out_edge = circ.get_nth_out_edge(v, 1);
      if (!circ.detect_final_Op(circ.target(c_out_edge)) ||
          circ.n_out_edges_of_type(v, EdgeType::Boolean) != 0) {
        if (dry_run) return {modified, false};
        throw CircuitInvalidity(
            "Cannot commute Measure through classical operations to the end "
            "of the circuit");
      }
      Edge out_edge = circ.get_nth_out_edge(v, 0);
      Edge current_edge = follow_until_noncommuting(circ, out_edge, true);
      Vertex current_vertex = circ.target(current_edge);
      // If we haven't moved it to an output, we can't continue
      if (!is_final_q_type(circ.get_OpType_from_Vertex(current_vertex))) {
        if (dry_run) return {modified, false};
        throw CircuitInvalidity(
            "Cannot commute Measure through quantum gates to the end of the "
            "circuit");
      }
      if (dry_run) continue;
      // If the measure was already at an output, nothing to do
      if (current_edge == out_edge) continue;
      Edge in_edge = circ.get_nth_in_edge(v, 0);
      // Rewire measure
      circ.add_edge(
          {circ.source(in_edge), circ.get_source_port(in_edge)},
          {circ.target(out_edge), circ.get_target_port(out_edge)},
          EdgeType::Quantum);
      circ.remove_edge(in_edge);
      circ.remove_edge(out_edge);
      circ.add_edge(
          {circ.source(current_edge), circ.get_source_port(current_edge)},
          {v, 0}, EdgeType::Quantum);
      circ.add_edge({v, 0}, {current_vertex, 0}, EdgeType::Quantum);
      circ.remove_edge(current_edge);
      modified = true;
    } else {
      Command com = command_from_vertex(circ, v);
      unit_set_t mesured_units;
      // Raise an error if there are any boxes or conditionals with internal
      // measures.
      if (check_only_end_measures(com, mesured_units)) {
        if (dry_run) return {modified, false};
        if (optype == OpType::Conditional) {
          throw CircuitInvalidity("Cannot delay measures inside a conditional");
        } else {
          throw CircuitInvalidity("Cannot delay measures inside a circuit box");
        }
      }
      // If the command has internal end measurements we allow them if they are
      // at the end of the full circuit. We can't delay the measures inside
      // boxes and conditionals.
      for (const auto& unit : mesured_units) {
        auto args = com.get_args();
        port_t out_port = std::distance(
            args.begin(), std::find(args.begin(), args.end(), unit));
        Edge current_edge = circ.get_nth_out_edge(v, out_port);
        // Follow the edge through swaps and boxes, but not through-commuting
        // gates
        current_edge = follow_until_noncommuting(circ, current_edge, false);
        Vertex current_vertex = circ.target(current_edge);
        if (!is_final_q_type(circ.get_OpType_from_Vertex(current_vertex))) {
          if (dry_run) return {modified, false};
          throw CircuitInvalidity(
              "Cannot commute Measure through quantum gates to the end of the "
              "circuit");
        }
      }
    }
  }
  return {modified, true};
}

}  // namespace DelayMeasures

}  // namespace Transforms

}  // namespace tket

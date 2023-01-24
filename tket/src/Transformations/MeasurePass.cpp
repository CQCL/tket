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

#include <tuple>

#include "Circuit/DAGDefs.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Transform.hpp"

namespace tket {

namespace Transforms {

Transform delay_measures() {
  return Transform(
      [](Circuit& circ) { return run_delay_measures_(circ, false).first; });
}

/** Commute all measurement gates to the end of the circuit.
 * @param circ The circuit to delay measurements in.
 * @param dry_run If true, do not modify the circuit, just check if it is
 * possible to delay.
 *
 * @throws CircuitInvalidity if it is not possible to delay and dry_run is
 * false.
 *
 * @return A pair of booleans. The first indicates whether the circuit was
 * changed, and the second indicates whether it was possible to delay (i.e.
 * there where no errors).
 **/
std::pair<bool, bool> run_delay_measures_(Circuit& circ, bool dry_run) {
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
      Edge current_edge = out_edge;
      Vertex current_vertex = circ.target(current_edge);
      port_t current_port = circ.get_target_port(current_edge);
      OpType current_optype = circ.get_OpType_from_Vertex(current_vertex);
      while (!is_final_q_type(current_optype) &&
             (current_optype == OpType::SWAP ||
              circ.commutes_with_basis(
                  current_vertex, Pauli::Z, PortType::Target, current_port))) {
        if (current_optype == OpType::SWAP) {
          // Update from SWAP
          current_edge =
              circ.get_nth_out_edge(current_vertex, 1 - current_port);
        } else {
          // Update to successor
          current_edge = circ.get_nth_out_edge(current_vertex, current_port);
        }
        current_vertex = circ.target(current_edge);
        current_port = circ.get_target_port(current_edge);
        current_optype = circ.get_OpType_from_Vertex(current_vertex);
        ;
      }
      // If we haven't moved it to an output, we can't continue
      if (!is_final_q_type(current_optype)) {
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
    } else if (optype == OpType::CircBox || optype == OpType::CustomGate) {
      // Raise an error if there are any boxes with internal measures.
      //
      // TODO We could (recursively) check if the measures are done at the
      // end of the internal circuits and not altered after that point, as
      // done in the NoMidMeasurePredicate. It should also check any
      // Conditionals.
      const Box& box = static_cast<const Box&>(*circ.get_Op_ptr_from_Vertex(v));
      for (const Command& c : *box.to_circuit()) {
        if (c.get_op_ptr()->get_type() == OpType::Measure) {
          if (dry_run) return {modified, false};
          throw CircuitInvalidity("Cannot delay measures inside a circuit box");
        }
      }
    }
  }
  return {modified, true};
}

}  // namespace Transforms

}  // namespace tket

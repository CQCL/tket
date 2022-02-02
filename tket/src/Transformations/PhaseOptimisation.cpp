// Copyright 2019-2022 Cambridge Quantum Computing
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

#include "PhaseOptimisation.hpp"

#include "Transform.hpp"

namespace tket {

namespace Transforms {

static void recursive_smash_CX_PhaseGadgets(
    Circuit &circ, Vertex &vert, VertexList &bin, bool &success);

Transform smash_CX_PhaseGadgets() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      recursive_smash_CX_PhaseGadgets(circ, v, bin, success);
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

static void recursive_smash_CX_PhaseGadgets(
    Circuit &circ, Vertex &vert, VertexList &bin, bool &success) {
  if (circ.get_OpType_from_Vertex(vert) == OpType::PhaseGadget) {
    for (unsigned e = 0; e < circ.n_in_edges(vert); e++) {
      Edge in_e = circ.get_nth_in_edge(vert, e);
      Edge out_e = circ.get_nth_out_edge(vert, e);
      Vertex prev_vert = circ.source(in_e);
      if (circ.get_OpType_from_Vertex(prev_vert) == OpType::CX &&
          circ.get_source_port(in_e) == 1) {
        Vertex next_vert = circ.target(out_e);
        if (circ.get_OpType_from_Vertex(next_vert) == OpType::CX &&
            circ.get_target_port(out_e) == 1) {
          Edge linker = circ.get_nth_in_edge(next_vert, 0);
          if (linker == circ.get_nth_out_edge(prev_vert, 0)) {
            success = true;
            circ.remove_edge(linker);
            port_t new_port = circ.n_in_edges(vert);
            circ.add_edge({prev_vert, 0}, {vert, new_port}, EdgeType::Quantum);
            circ.add_edge({vert, new_port}, {next_vert, 0}, EdgeType::Quantum);
            VertexList to_detach = {prev_vert, next_vert};
            bin.push_back(prev_vert);
            bin.push_back(next_vert);
            circ.remove_vertices(
                to_detach, Circuit::GraphRewiring::Yes,
                Circuit::VertexDeletion::No);
            e--;  // continue to add any more from the current port
          }
        }
      }
    }
    // Fix up the qubit count for the phase gadget we inserted.
    std::vector<Expr> v_params =
        circ.get_Op_ptr_from_Vertex(vert)->get_params();
    circ.dag[vert].op =
        get_op_ptr(OpType::PhaseGadget, v_params, circ.n_in_edges(vert));
  }
}

static bool align_phases_all(Circuit &circ) {
  bool success = false;
  VertexList bin;
  VertexVec vertices = circ.vertices_in_order();
  auto index = boost::get(boost::vertex_index, circ.dag);
  for (const Vertex &v : vertices) {
    if (circ.get_OpType_from_Vertex(v) == OpType::PhaseGadget) {
      std::map<Vertex, std::map<port_t, port_t>>
          parent_to_port_map;  // maps ports on initial phase to ports on each
                               // parent
      EdgeVec i_edges = circ.get_in_edges(v);
      unsigned n_qubits = i_edges.size();
      for (EdgeVec::iterator e = i_edges.begin(); e != i_edges.end(); ++e) {
        Edge ed = *e;
        Vertex source = boost::source(ed, circ.dag);
        OpType type = circ.get_OpType_from_Vertex(source);
        while (!(type == OpType::PhaseGadget || is_initial_q_type(type))) {
          ed = circ.get_last_edge(source, ed);
          source = boost::source(ed, circ.dag);
          type = circ.get_OpType_from_Vertex(source);
        }
        if (is_initial_q_type(type) || circ.get_source_port(ed) > n_qubits) {
          continue;
        }
        std::map<Vertex, std::map<port_t, port_t>>::iterator par =
            parent_to_port_map.find(source);
        if (par == parent_to_port_map.end()) {
          std::map<port_t, port_t> port_map;
          port_map.insert({circ.get_target_port(*e), circ.get_source_port(ed)});
          parent_to_port_map.insert({source, port_map});
        } else {
          par->second.insert(
              {circ.get_target_port(*e), circ.get_source_port(ed)});
        }
      }
      // build replacement as disconnected graph
      Circuit phase_replacement;
      std::map<port_t, Vertex> port_to_input;
      std::map<port_t, Vertex> port_to_output;
      std::map<port_t, bool> port_remaining;
      for (port_t p = 0; p < n_qubits; p++) {
        Vertex in = phase_replacement.add_vertex(OpType::Input);
        Vertex out = phase_replacement.add_vertex(OpType::Output);
        phase_replacement.boundary.insert({Qubit(p), in, out});
        port_to_input.insert({p, in});
        port_to_output.insert({p, out});
        port_remaining.insert({p, true});
      }
      Vertex gadget = phase_replacement.add_vertex(
          phase_replacement.get_Op_ptr_from_Vertex(v));
      // collect incoming edges by parent and add where preference is possible
      typedef std::function<bool(
          std::pair<Vertex, std::map<port_t, port_t>>,
          std::pair<Vertex, std::map<port_t, port_t>>)>
          Comp;
      std::set<std::pair<Vertex, std::map<port_t, port_t>>, Comp> sorted_p2p(
          parent_to_port_map.begin(), parent_to_port_map.end(),
          [&](std::pair<Vertex, std::map<port_t, port_t>> const &a,
              std::pair<Vertex, std::map<port_t, port_t>> const &b) {
            if (a.second.size() == b.second.size()) {
              return index[a.first] < index[b.first];
            }
            return a.second.size() > b.second.size();
          });
      for (std::pair<Vertex, std::map<port_t, port_t>> par : sorted_p2p) {
        for (std::map<port_t, port_t>::iterator port_pair = par.second.begin();
             port_pair != par.second.end(); ++port_pair) {
          if (port_remaining[port_pair->second]) {
            phase_replacement.add_edge(
                {port_to_input[port_pair->first], 0},
                {gadget, port_pair->second}, EdgeType::Quantum);
            phase_replacement.add_edge(
                {gadget, port_pair->second},
                {port_to_output[port_pair->first], 0}, EdgeType::Quantum);
            port_remaining[port_pair->second] = false;
            port_to_input.erase(port_pair->first);
          }
        }
      }
      // add edges for remaining ports
      port_t spare_port = 0;
      for (std::map<port_t, Vertex>::iterator in = port_to_input.begin();
           in != port_to_input.end(); ++in) {
        while (!port_remaining[spare_port]) spare_port++;
        phase_replacement.add_edge(
            {in->second, 0}, {gadget, spare_port}, EdgeType::Quantum);
        phase_replacement.add_edge(
            {gadget, spare_port}, {port_to_output[in->first], 0},
            EdgeType::Quantum);
        spare_port++;
      }
      // substitute into circuit
      Subcircuit sub = {
          {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
      bin.push_back(v);
      circ.substitute(phase_replacement, sub, Circuit::VertexDeletion::No);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

Transform align_PhaseGadgets() { return Transform(align_phases_all); }

}  // namespace Transforms

}  // namespace tket

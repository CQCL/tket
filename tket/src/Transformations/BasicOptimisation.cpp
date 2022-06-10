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

#include "BasicOptimisation.hpp"

#include <optional>

#include "Characterisation/DeviceCharacterisation.hpp"
#include "Characterisation/ErrorTypes.hpp"
#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Command.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Decomposition.hpp"
#include "Gate/Gate.hpp"
#include "Gate/GatePtr.hpp"
#include "Gate/Rotation.hpp"
#include "Transform.hpp"
#include "Utils/Assert.hpp"
#include "Utils/EigenConfig.hpp"

namespace tket {

namespace Transforms {

static bool redundancy_removal(Circuit &circ);
static bool remove_redundancy(
    Circuit &circ, const Vertex &vert, VertexList &bin,
    std::set<IVertex> &new_affected_verts, IndexMap &im);
static bool commute_singles_to_front(Circuit &circ);

Transform remove_redundancies() { return Transform(redundancy_removal); }

// this method annihilates all primitives next to each other (accounting for
// previous annihilations)
// also removes redundant non-classically controlled Z basis gates before a z
// basis measurement so that eg. -H-X-X-H- always annihilates to -----
static bool redundancy_removal(Circuit &circ) {
  bool success = false;
  bool found_redundancy = true;
  IndexMap im = circ.index_map();
  std::set<IVertex> old_affected_verts;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    old_affected_verts.insert({im.at(v), v});
  }
  VertexList bin;
  while (found_redundancy) {
    std::set<IVertex> new_affected_verts;
    for (const IVertex &p : old_affected_verts) {
      remove_redundancy(circ, p.second, bin, new_affected_verts, im);
    }
    found_redundancy = new_affected_verts.size() != 0;
    success |= found_redundancy;
    old_affected_verts = new_affected_verts;
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

// called by the previous method. This should generally not be called
// independently
static bool remove_redundancy(
    Circuit &circ, const Vertex &vert, VertexList &bin,
    std::set<IVertex> &new_affected_verts, IndexMap &im) {
  const Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
  const OpDesc desc = op->get_desc();
  if (!desc.is_gate()) return false;
  if (circ.n_out_edges(vert) == 0 || circ.n_in_edges(vert) == 0) {
    return false;  // either a boundary vert or we have already removed it
  }

  auto remove_single_vertex = [&bin, &circ, &new_affected_verts,
                               &im](const Vertex &v_remove) {
    bin.push_back(v_remove);
    for (const Vertex &l : circ.get_predecessors(v_remove)) {
      new_affected_verts.insert({im.at(l), l});
    }
    circ.remove_vertex(
        v_remove, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
  };
  // remove 0 angle rotations from circuit
  std::optional<double> a = op->is_identity();
  if (a) {
    remove_single_vertex(vert);
    circ.add_phase(a.value());
    return true;
  } else if (desc.type() == OpType::noop) {
    // remove "noop" gates from circuit
    remove_single_vertex(vert);
    return true;
  }
  VertexVec kids = circ.get_successors(vert);

  // if op is immediately followed by a z basis measurement on all qubits,
  // remove
  if (circ.n_out_edges_of_type(vert, EdgeType::Classical) == 0) {
    bool z_followed_by_measures = true;
    for (port_t port = 0; port < kids.size() && z_followed_by_measures;
         port++) {
      if (circ.get_OpType_from_Vertex(kids[port]) == OpType::Measure) {
        z_followed_by_measures &=
            circ.commutes_with_basis(vert, Pauli::Z, PortType::Source, port);
      } else {
        z_followed_by_measures = false;
      }
    }
    if (z_followed_by_measures) {
      remove_single_vertex(vert);
      return true;
    }
  }

  // check that both the vertex and its successor have each other and only each
  // other
  if ((kids.size() == 1) && (circ.get_predecessors(kids[0]).size() == 1)) {
    // check that the ports match up between vertices
    Vertex b = kids[0];
    EdgeVec ins = circ.get_in_edges(b);
    for (const Edge &in : ins) {
      if (circ.get_source_port(in) != circ.get_target_port(in)) return false;
    }

    // check that the classical edges match up correctly
    if (circ.n_in_edges_of_type(vert, EdgeType::Boolean) != 0) return false;

    const Op_ptr b_op = circ.get_Op_ptr_from_Vertex(b);
    const OpDesc b_desc = b_op->get_desc();

    if (!b_desc.is_oneway()) {
      // if A = B.dagger(), AB = I
      // This method cannot detect matches between rotation gates.
      // Rotation gates are covered by the rotation gate combiner, everything
      // else in this.
      if (*b_op->dagger() == *op) {
        bin.push_back(vert);
        bin.push_back(b);
        VertexVec last_verts = circ.get_predecessors(vert);
        for (const Vertex &l : last_verts) {
          new_affected_verts.insert({im.at(l), l});
        }
        VertexList to_detach{vert, b};
        // detached from circuit but not removed from graph
        circ.remove_vertices(
            to_detach, Circuit::GraphRewiring::Yes,
            Circuit::VertexDeletion::No);
        return true;
      } else if (desc.is_rotation()) {
        // combine two rotation gates together, then if the combined
        // operation is the identity up to phase, remove from circuit
        if (b_desc.type() == desc.type()) {
          Expr expr1 = op->get_params()[0];
          Expr expr2 = b_op->get_params()[0];
          VertexVec last_verts = circ.get_predecessors(vert);
          for (const Vertex &l : last_verts) {
            new_affected_verts.insert({im.at(l), l});
          }
          circ.remove_vertex(
              b, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
          bin.push_back(b);
          std::vector<Expr> params_new = {expr1 + expr2};
          Op_ptr op_new = get_op_ptr(desc.type(), params_new, ins.size());
          std::optional<double> a = op_new->is_identity();
          if (a) {
            bin.push_back(vert);
            circ.remove_vertex(
                vert, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            circ.add_phase(a.value());
          } else {
            new_affected_verts.insert({im[vert], vert});
            circ.dag[vert].op = op_new;
          }
          return true;
        }
      }
    }
  }
  return false;
}

Transform commute_through_multis() {
  return Transform(commute_singles_to_front);
}

// moves single qubit operations past multiqubit operations they commute with,
// towards front of circuit (hardcoded)
static bool commute_singles_to_front(Circuit &circ) {
  bool success = false;
  // follow each qubit path from output to input
  for (const Qubit &q : circ.all_qubits()) {
    Vertex prev_v = circ.get_out(q);
    Edge current_e = circ.get_nth_in_edge(prev_v, 0);
    Vertex current_v = circ.source(current_e);
    while (!is_initial_q_type(circ.get_OpType_from_Vertex(current_v))) {
      const Op_ptr curr_op = circ.get_Op_ptr_from_Vertex(current_v);
      // if current vertex is a multiqubit gate
      if (circ.n_in_edges_of_type(current_v, EdgeType::Quantum) > 1 &&
          curr_op->get_desc().is_gate()) {
        // need gate check to access commuting_basis
        const std::pair<port_t, port_t> ports = circ.get_ports(current_e);
      check_prev_commutes:
        const Op_ptr prev_op = circ.get_Op_ptr_from_Vertex(prev_v);
        // if previous vertex is single qubit gate
        if (prev_op->get_desc().is_gate() &&
            circ.n_in_edges_of_type(prev_v, EdgeType::Quantum) == 1) {
          const std::optional<Pauli> prev_colour =
              circ.commuting_basis(prev_v, PortType::Target, ports.second);

          if (circ.commutes_with_basis(
                  current_v, prev_colour, PortType::Source, ports.first)) {
            // subsequent op on qubit path is a single qubit gate
            // and commutes with current multi qubit gate
            success = true;
            circ.remove_vertex(
                prev_v, Circuit::GraphRewiring::Yes,
                Circuit::VertexDeletion::No);
            Edge rewire_edge = circ.get_nth_in_edge(current_v, ports.first);
            circ.rewire(prev_v, {rewire_edge}, {EdgeType::Quantum});
            current_e = circ.get_nth_out_edge(current_v, ports.first);
            prev_v = circ.target(current_e);
            // check if new previous gate can be commuted through too
            goto check_prev_commutes;
          }
        }
      }
      // move to next vertex (towards input)
      prev_v = current_v;
      std::tie(current_v, current_e) = circ.get_prev_pair(current_v, current_e);
    }
  }

  return success;
}

// helper class subcircuits representing 2qb interactions
struct Interaction {
  Interaction(const Qubit &_q0, const Qubit &_q1) : q0(_q0), q1(_q1) {}
  Qubit q0;  // Qubit numbers
  Qubit q1;
  Edge e0;  // In edges starting interaction
  Edge e1;
  unsigned count;      // Number of two qubit gates in interaction
  VertexSet vertices;  // Vertices in interaction subcircuit
};

static bool replace_two_qubit_interaction(
    Circuit &circ, Interaction &i, std::map<Qubit, Edge> &current_edges,
    VertexList &bin, const double cx_fidelity = 1.) {
  EdgeVec in_edges = {i.e0, i.e1};
  EdgeVec out_edges = {current_edges[i.q0], current_edges[i.q1]};
  Edge next0, next1;
  bool q0_is_out = is_final_q_type(
      circ.get_OpType_from_Vertex(circ.target(current_edges[i.q0])));
  bool q1_is_out = is_final_q_type(
      circ.get_OpType_from_Vertex(circ.target(current_edges[i.q1])));
  if (!q0_is_out) {
    next0 = circ.get_next_edge(
        circ.target(current_edges[i.q0]), current_edges[i.q0]);
  }
  if (!q1_is_out) {
    next1 = circ.get_next_edge(
        circ.target(current_edges[i.q1]), current_edges[i.q1]);
  }
  Subcircuit sub = {in_edges, out_edges, i.vertices};
  Circuit subc = circ.subcircuit(sub);
  Eigen::Matrix4cd mat = get_matrix_from_2qb_circ(subc);
  Circuit replacement = two_qubit_canonical(mat);
  // TODO: for now we decompose all the way to CX. In the future, this pass
  // should output TK2, and decompose to CX (or other gates) later if necessary.
  TwoQbFidelities fid;
  fid.CX_fidelity = cx_fidelity;
  (decompose_TK2(fid) >> squash_1qb_to_tk1()).apply(replacement);
  const int nb_cx_old = subc.count_gates(OpType::CX);
  const int nb_cx_new = replacement.count_gates(OpType::CX);
  if (nb_cx_new < nb_cx_old) {
    bin.insert(bin.end(), sub.verts.begin(), sub.verts.end());
    circ.substitute(replacement, sub, Circuit::VertexDeletion::No);
    if (!q0_is_out) {
      current_edges[i.q0] = circ.get_last_edge(circ.source(next0), next0);
    }
    if (!q1_is_out) {
      current_edges[i.q1] = circ.get_last_edge(circ.source(next1), next1);
    }
    return true;
  } else {
    return false;
  }
}

Transform commute_and_combine_HQS2() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      EdgeVec outs = circ.get_all_out_edges(v);
      if (circ.get_OpType_from_Vertex(v) == OpType::ZZMax && outs.size() == 2) {
        Vertex next0 = boost::target(outs[0], circ.dag);
        Vertex next1 = boost::target(outs[1], circ.dag);
        if (next0 == next1 &&
            circ.get_OpType_from_Vertex(next0) == OpType::ZZMax) {
          success = true;
          EdgeVec h_in = circ.get_in_edges(v);
          EdgeVec h_out = circ.get_all_out_edges(next0);
          if (circ.get_target_port(outs[0]) != 0) {
            h_out = {h_out[1], h_out[0]};
          }
          bin.push_back(v);
          bin.push_back(next0);
          Subcircuit sub = {h_in, h_out};
          circ.substitute(
              CircPool::two_Rz1(), sub, Circuit::VertexDeletion::No);
          circ.add_phase(0.5);
          continue;
        }
        if (circ.get_OpType_from_Vertex(next0) == OpType::Rz) {
          success = true;
          circ.remove_vertex(
              next0, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
          Edge in_0 = circ.get_nth_in_edge(v, 0);
          circ.rewire(next0, {in_0}, {EdgeType::Quantum});
        }
        if (circ.get_OpType_from_Vertex(next1) == OpType::Rz) {
          success = true;
          circ.remove_vertex(
              next1, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
          Edge in_1 = circ.get_nth_in_edge(v, 1);
          circ.rewire(next1, {in_1}, {EdgeType::Quantum});
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

// TODO:: Work around classically controlled stuff
Transform two_qubit_squash(double cx_fidelity) {
  return Transform([cx_fidelity](Circuit &circ) {
    bool success = false;
    VertexList bin;
    // Get map from vertex/port to qubit number
    std::map<VertPort, Qubit> v_to_qb;
    std::map<Qubit, Edge> current_edge_on_qb;
    std::vector<Interaction> i_vec;
    std::map<Qubit, int> current_interaction;
    for (const Qubit &qb : circ.all_qubits()) {
      for (const VertPort &vp : circ.unit_path(qb)) {
        v_to_qb.insert({vp, qb});
      }
      Vertex input = circ.get_in(qb);
      Edge e = circ.get_nth_out_edge(input, 0);
      current_edge_on_qb[qb] = e;
      current_interaction[qb] = -1;
    }
    SliceVec slices = circ.get_slices();
    slices.insert(slices.begin(), circ.q_inputs());
    slices.push_back(circ.q_outputs());
    for (SliceVec::iterator s = slices.begin(); s != slices.end(); ++s) {
      for (Slice::iterator v = s->begin(); v != s->end(); ++v) {
        const Op_ptr o = circ.get_Op_ptr_from_Vertex(*v);
        OpType type = o->get_type();
        unsigned n_ins = circ.n_in_edges_of_type(*v, EdgeType::Quantum);
        // Measures, resets, outputs, barriers, symbolic gates, and many-qubit
        // gates close interactions
        if (is_projective_type(type) || is_final_q_type(type) ||
            type == OpType::Barrier || n_ins > 2 ||
            !o->free_symbols().empty()) {
          for (port_t port = 0; port < n_ins; port++) {
            Qubit q = v_to_qb.at({*v, port});
            int i = current_interaction[q];
            if (i != -1) {
              if (i_vec[i].count >= 2) {
                // Replace subcircuit
                success |= replace_two_qubit_interaction(
                    circ, i_vec[i], current_edge_on_qb, bin, cx_fidelity);
              }
              current_interaction[i_vec[i].q0] = -1;
              current_interaction[i_vec[i].q1] = -1;
            }
            if (!is_final_q_type(type)) {
              current_edge_on_qb[q] =
                  circ.get_next_edge(*v, current_edge_on_qb[q]);
            }
          }
          continue;
        }

        // Check for 2qb gate
        if (circ.n_in_edges_of_type(*v, EdgeType::Quantum) == 2) {
          Qubit q0 = v_to_qb.at({*v, 0});
          Qubit q1 = v_to_qb.at({*v, 1});
          int i0 = current_interaction[q0];
          int i1 = current_interaction[q1];
          // If they are already interacting, extend it
          if (i0 != -1 && i0 == i1) {
            i_vec[i0].count++;
            i_vec[i0].vertices.insert(*v);
            current_edge_on_qb[q0] =
                circ.get_next_edge(*v, current_edge_on_qb[q0]);
            current_edge_on_qb[q1] =
                circ.get_next_edge(*v, current_edge_on_qb[q1]);
            continue;
          }
          // End any other interactions on q0
          if (i0 != -1) {
            if (i_vec[i0].count >= 2) {
              // Replace subcircuit
              success |= replace_two_qubit_interaction(
                  circ, i_vec[i0], current_edge_on_qb, bin, cx_fidelity);
            }
            current_interaction[i_vec[i0].q0] = -1;
            current_interaction[i_vec[i0].q1] = -1;
          }
          // End any other interactions on q1
          if (i1 != -1) {
            if (i_vec[i1].count >= 2) {
              // Replace subcircuit

              success |= replace_two_qubit_interaction(
                  circ, i_vec[i1], current_edge_on_qb, bin, cx_fidelity);
            }
            current_interaction[i_vec[i1].q0] = -1;
            current_interaction[i_vec[i1].q1] = -1;
          }
          // Add new interaction
          Interaction new_i(q0, q1);
          new_i.e0 = current_edge_on_qb[q0];
          new_i.e1 = current_edge_on_qb[q1];
          new_i.count = 1;
          new_i.vertices = {*v};
          current_interaction[q0] = i_vec.size();
          current_interaction[q1] = i_vec.size();
          i_vec.push_back(new_i);
          current_edge_on_qb[q0] =
              circ.get_next_edge(*v, current_edge_on_qb[q0]);
          current_edge_on_qb[q1] =
              circ.get_next_edge(*v, current_edge_on_qb[q1]);
          continue;
        }

        // Otherwise, we don't care about other vertices, so just update edges
        // and add vertices if interactions exist
        for (port_t i = 0; i < circ.n_in_edges(*v); i++) {
          Qubit q = v_to_qb.at({*v, i});
          current_edge_on_qb[q] = circ.get_next_edge(*v, current_edge_on_qb[q]);
          int inter = current_interaction[q];
          if (inter != -1) {
            i_vec[inter].vertices.insert(*v);
          }
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

// Given a 'SWAP_chain', finds edge in chain (or qubit wire) with best fidelity
// and rewires the associated single qubit vertex in to it
static bool find_edge_rewire_vertex(
    Circuit &circ,
    std::pair<std::vector<std::pair<Edge, double>>, Vertex> &entry) {
  std::vector<std::pair<Edge, double>> candidates = entry.first;
  std::pair<Edge, double> best_pair = candidates[0];
  bool do_rewire = false;  // prevents rewiring if current edge is best edge
                           // (causes many issues...)
  for (unsigned i = 1; i < candidates.size(); i++) {  // find best edge
    if (candidates[i].second > best_pair.second) {
      best_pair = candidates[i];
      do_rewire = true;
    }
  }
  if (do_rewire) {
    // rewire in to best edge
    circ.remove_vertex(
        entry.second, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
    circ.rewire(entry.second, {best_pair.first}, {EdgeType::Quantum});
  }
  return do_rewire;
}

// Given a SWAP vertex, which has some predecessor SWAP vertex, find the
// 'SWAP_chain' this predecessor SWAP vertex is in and add to it
static void extend_SWAP_chain(
    std::list<std::pair<std::vector<std::pair<Edge, double>>, Vertex>>
        &swap_chains,
    Edge entry_edge, Node entry_node, const Edge &match, const Circuit &circ,
    const DeviceCharacterisation &characterisation) {
  for (auto it = swap_chains.begin(); it != swap_chains.end(); ++it) {
    if ((*it).first.back().first == match) {
      // extend chain, adding a new pair of edge and double to the end
      (*it).first.push_back(
          {entry_edge,
           1.0 - characterisation.get_error(
                     entry_node, circ.get_OpType_from_Vertex((*it).second))});
      return;
    }
  }
}

// Finds sequences of adjacent SWAP gates, with a predecessor Single Qubit
// Vertex. The error rate of the required Single Qubit Vertex is stored for each
// of the PhysicalQubits the LogicalQubit passes through Once 'SWAP_chains' are
// found throughout the whole circuit, predecessor Single Qubit Vertices are
// rewired into the edge with best error rate
static bool find_rewire_sq(
    Circuit &circ, const DeviceCharacterisation &characterisation) {
  std::list<std::pair<std::vector<std::pair<Edge, double>>, Vertex>>
      swap_chains;
  for (auto it = circ.begin(); it != circ.end(); ++it) {
    if (it->get_op_ptr()->get_type() == OpType::SWAP) {
      // find SWAP, if either predecessor is a single qubit unitary
      // find resulting swap chain...
      Vertex swap_vert = it.get_vertex();
      unit_vector_t qubits = it->get_args();
      node_vector_t nodes = {qubits.begin(), qubits.end()};
      VertexVec pred_verts = circ.get_predecessors(swap_vert);
      EdgeVec pred_edges = circ.get_in_edges(swap_vert);
      EdgeVec post_edges = circ.get_all_out_edges(swap_vert);
      for (unsigned i = 0; i < pred_verts.size(); i++) {
        OpType optype = circ.get_OpType_from_Vertex(
            pred_verts[i]);  // OPT IS PREDECESSOR OPTYPE
        if (circ.detect_singleq_unitary_op(pred_verts[i])) {
          // wire has single qubit unitary->add a new swap chain
          std::vector<std::pair<Edge, double>> swap_chain = {
              {pred_edges[i],
               1.0 - characterisation.get_error(nodes[i], optype)},
              {post_edges[1 - i],
               1.0 - characterisation.get_error(nodes[1 - i], optype)}};
          swap_chains.push_back(std::make_pair(swap_chain, pred_verts[i]));
        } else if (optype == OpType::SWAP) {  // wire has swap->assume this swap
                                              // already in chain, find chain
          Edge pre_edge = pred_edges[i];
          extend_SWAP_chain(
              swap_chains, post_edges[1 - i], nodes[1 - i], pre_edge, circ,
              characterisation);
        }
      }
    }
  }
  // having produced swap chains, now find best qubit for gate to act on and
  // implement
  bool success = false;
  while (!swap_chains.empty()) {
    if (find_edge_rewire_vertex(circ, (*swap_chains.begin())) == true) {
      success = true;
    }
    swap_chains.erase(swap_chains.begin());
  }
  return success;
}

static Transform commute_SQ_gates_through_SWAPS_helper(
    const DeviceCharacterisation &characterisation) {
  return Transform([characterisation](Circuit &circ) {
    bool success = false;
    while (find_rewire_sq(circ, characterisation)) {
      success = true;
    }
    return success;
  });
}
Transform commute_SQ_gates_through_SWAPS(const avg_node_errors_t &node_errors) {
  return commute_SQ_gates_through_SWAPS_helper(
      DeviceCharacterisation(node_errors));
}
Transform commute_SQ_gates_through_SWAPS(const op_node_errors_t &node_errors) {
  return commute_SQ_gates_through_SWAPS_helper(
      DeviceCharacterisation(node_errors));
}

Transform absorb_Rz_NPhasedX() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexSet all_bins;

    // Start by squashing Rz gates
    success |= squash_1qb_to_pqp(OpType::Rz, OpType::Rx).apply(circ);

    // Loop through all NPhasedX gates
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      if (op->get_type() == OpType::NPhasedX) {
        // gather surrounding Rz gates
        unsigned arity = op->n_qubits();
        std::vector<Expr> in_rz(arity);
        std::vector<Expr> out_rz(arity);
        EdgeVec in_edges = circ.get_in_edges_of_type(v, EdgeType::Quantum);
        EdgeVec out_edges = circ.get_out_edges_of_type(v, EdgeType::Quantum);
        TKET_ASSERT(in_edges.size() == arity);
        TKET_ASSERT(out_edges.size() == arity);
        for (unsigned i = 0; i < arity; ++i) {
          Vertex in_v = circ.source(in_edges[i]);
          Op_ptr in_op = circ.get_Op_ptr_from_Vertex(in_v);
          Vertex out_v = circ.target(out_edges[i]);
          Op_ptr out_op = circ.get_Op_ptr_from_Vertex(out_v);

          if (in_op->get_type() == OpType::Rz) {
            in_rz[i] = -in_op->get_params().at(0);
          } else {
            in_rz[i] = 0.;
          }
          if (out_op->get_type() == OpType::Rz) {
            out_rz[i] = out_op->get_params().at(0);
          } else {
            out_rz[i] = 0.;
          }
        }

        // Find out which Rz angle is most popular.
        // Note that we only compare expr[i] with expr[j] when j < i. This means
        // that only the largest i from a set of equivalent exprs will have the
        // right occurence count, but that is good enough.
        std::vector<Expr> all_rz = in_rz;
        all_rz.insert(all_rz.end(), out_rz.begin(), out_rz.end());
        std::vector<unsigned> occurences_count(2 * arity);
        for (unsigned i = 0; i < 2 * arity; ++i) {
          unsigned cnt = 0;
          for (unsigned j = 0; j < i; ++j) {
            if (equiv_expr(all_rz[i], all_rz[j], 4)) {
              ++cnt;
            }
          }
          occurences_count[i] = cnt;
        }
        unsigned max_i =
            std::max_element(occurences_count.begin(), occurences_count.end()) -
            occurences_count.begin();
        Expr absorb_rz = all_rz[max_i];

        if (!equiv_0(absorb_rz, 4)) {
          success = true;

          // Subtract absorb_rz in NPhasedX
          std::vector<Expr> new_params = op->get_params();
          TKET_ASSERT(new_params.size() == 2);
          new_params[1] += absorb_rz;
          circ.dag[v] = get_op_ptr(OpType::NPhasedX, new_params, arity);

          // Finally, adjust +-absorb_rz in Rz everywhere around
          for (unsigned i = 0; i < arity; ++i) {
            Vertex in_v = circ.source(in_edges[i]);
            Op_ptr in_op = circ.get_Op_ptr_from_Vertex(in_v);
            Vertex out_v = circ.target(out_edges[i]);
            Op_ptr out_op = circ.get_Op_ptr_from_Vertex(out_v);

            Expr angle;
            Edge in_e, out_e;
            VertexSet bin;
            if (in_op->get_type() == OpType::Rz) {
              angle = in_op->get_params().at(0) + absorb_rz;
              out_e = in_edges[i];
              in_e = circ.get_last_edge(in_v, out_e);
              bin = {in_v};
            } else {
              angle = absorb_rz;
              out_e = in_edges[i];
              in_e = out_e;
              bin = {};
            }
            Subcircuit sub{{in_e}, {out_e}, bin};
            Circuit c(1);
            if (!equiv_0(angle, 4)) {
              c.add_op<unsigned>(OpType::Rz, angle, {0});
            }
            circ.substitute(c, sub, Circuit::VertexDeletion::No);
            all_bins.insert(bin.begin(), bin.end());

            if (out_op->get_type() == OpType::Rz) {
              angle = out_op->get_params().at(0) - absorb_rz;
              in_e = out_edges[i];
              out_e = circ.get_next_edge(out_v, in_e);
              bin = {out_v};
            } else {
              angle = -absorb_rz;
              in_e = out_edges[i];
              out_e = in_e;
              bin = {};
            }
            sub = Subcircuit{{in_e}, {out_e}, bin};
            c = Circuit(1);
            if (!equiv_0(angle, 4)) {
              c.add_op<unsigned>(OpType::Rz, angle, {0});
            }
            circ.substitute(c, sub, Circuit::VertexDeletion::No);
            all_bins.insert(bin.begin(), bin.end());
          }
        }
      }
    }
    circ.remove_vertices(
        all_bins, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    return success;
  });
}

}  // namespace Transforms

}  // namespace tket

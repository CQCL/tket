// Copyright Quantinuum
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

#include "tket/Transformations/BasicOptimisation.hpp"

#include <cmath>
#include <optional>
#include <tkassert/Assert.hpp>
#include <vector>

#include "tket/Characterisation/DeviceCharacterisation.hpp"
#include "tket/Characterisation/ErrorTypes.hpp"
#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/Transform.hpp"
#include "tket/Utils/Expression.hpp"
#include "tket/Utils/MatrixAnalysis.hpp"
#include "tket/Utils/UnitID.hpp"

namespace tket::Transforms {

static bool commute_singles_to_front(Circuit &circ);

// called by the previous method. This should generally not be called
// independently

Transform commute_through_multis() {
  return Transform(commute_singles_to_front);
}

// whether source and target commute
static bool ends_commute(const Circuit &circ, const Edge &e) {
  const std::pair<port_t, port_t> ports = circ.get_ports(e);
  const Vertex source = circ.source(e);
  const Vertex target = circ.target(e);

  // We currently do not support commuting multi-qubit gates
  // TODO: It would be useful to support commuting single-qubit gates with
  // classical conditioning
  if (circ.n_in_edges(source) > 1 && circ.n_in_edges(target) > 1) {
    return false;
  }

  auto colour = circ.commuting_basis(target, ports.second);
  return circ.commutes_with_basis(source, colour, ports.first);
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
      if (circ.n_in_edges_of_type(current_v, EdgeType::Quantum) > 1) {
        while (circ.n_in_edges_of_type(prev_v, EdgeType::Quantum) == 1 &&
               ends_commute(circ, current_e)) {
          // subsequent op on qubit path is a single qubit gate
          // and commutes with current multi qubit gate
          success = true;
          EdgeVec rewire_edges;
          op_signature_t edge_types;
          for (const Edge &e : circ.get_in_edges(prev_v)) {
            EdgeType type = circ.get_edgetype(e);
            Edge boundary_edge;
            // Currently, only purely-quantum operations can be commuted
            // through. This is guaranteed by `ends_commute`. It follows that
            // any wire out of `prev_v` must be EdgeType::Quantum.
            TKET_ASSERT(type == EdgeType::Quantum);
            boundary_edge = circ.get_last_edge(current_v, current_e);
            rewire_edges.push_back(boundary_edge);
            edge_types.push_back(type);
          }
          const port_t backup_port = circ.get_source_port(current_e);
          circ.remove_vertex(
              prev_v, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
          circ.rewire(prev_v, rewire_edges, edge_types);
          current_e = circ.get_nth_out_edge(current_v, backup_port);
          prev_v = circ.target(current_e);
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
    VertexList &bin, OpType target, double cx_fidelity, bool allow_swaps) {
  EdgeVec in_edges = {i.e0, i.e1};
  std::vector<std::optional<Edge>> out_edges = {
      current_edges[i.q0], current_edges[i.q1]};
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
  // Circuit to (potentially) substitute
  Subcircuit sub = {in_edges, out_edges, {}, i.vertices};
  Circuit subc = circ.subcircuit(sub);

  // Try to simplify using KAK
  Circuit replacement = subc;
  decompose_multi_qubits_TK2().apply(replacement);
  Eigen::Matrix4cd mat = get_matrix_from_2qb_circ(replacement);
  replacement = two_qubit_canonical(mat);
  TwoQbFidelities fid;
  fid.CX_fidelity = cx_fidelity;
  if (target != OpType::TK2) {
    decompose_TK2(fid, allow_swaps).apply(replacement);
  }

  // Whether to substitute old circuit with new
  bool substitute = false;
  for (Vertex v : subc.vertices_in_order()) {
    unsigned n_ins = subc.n_in_edges_of_type(v, EdgeType::Quantum);
    if (n_ins == 2 && subc.get_OpType_from_Vertex(v) != target) {
      // Old circuit has non-target gates => we need to substitute
      substitute = true;
      break;
    }
  }
  if (!substitute) {
    if (target == OpType::CX) {
      unsigned nb_2qb_old = subc.count_gates(target);
      unsigned nb_2qb_new = replacement.count_gates(target);
      substitute |= nb_2qb_new < nb_2qb_old;
    } else if (target == OpType::TK2) {
      unsigned cnt_2qb = 0;
      for (Vertex v : subc.vertices_in_order()) {
        unsigned n_ins = subc.n_in_edges_of_type(v, EdgeType::Quantum);
        cnt_2qb += n_ins == 2;
      }
      substitute |= cnt_2qb >= 2;
    }
  }

  if (substitute) {
    // Substitute interaction with new circuit
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
    // Leave circuit untouched
    return false;
  }
}

Transform two_qubit_squash(bool allow_swaps) {
  return two_qubit_squash(OpType::CX, 1., allow_swaps);
}

Transform two_qubit_squash(
    OpType target_2qb_gate, double cx_fidelity, bool allow_swaps) {
  const std::set<OpType> accepted_ots{OpType::CX, OpType::TK2};
  if (!accepted_ots.contains(target_2qb_gate)) {
    throw BadOpType(
        "KAKDecomposition currently supports CX and TK2. "
        "Cannot decompose to",
        target_2qb_gate);
  }
  if (cx_fidelity < 0 || cx_fidelity > 1) {
    throw std::invalid_argument("The CX fidelity must be between 0 and 1.");
  }

  return Transform([target_2qb_gate, cx_fidelity, allow_swaps](Circuit &circ) {
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
        // Ignore classical ops
        if (is_classical_type(type)) {
          continue;
        } else if (
            is_projective_type(type) || is_final_q_type(type) ||
            type == OpType::Barrier || type == OpType::Conditional ||
            n_ins > 2 || !o->free_symbols().empty()) {
          // Measures, resets, outputs, barriers, symbolic gates, conditionals
          // and many-qubit gates close interactions
          EdgeVec q_edges = circ.get_in_edges_of_type(*v, EdgeType::Quantum);
          std::vector<port_t> q_ports;
          for (const Edge &e : q_edges) {
            q_ports.push_back(circ.get_target_port(e));
          }
          for (const port_t &port : q_ports) {
            Qubit q = v_to_qb.at({*v, port});
            int i = current_interaction[q];
            if (i != -1) {
              if (i_vec[i].count >= 2) {
                // Replace subcircuit
                success |= replace_two_qubit_interaction(
                    circ, i_vec[i], current_edge_on_qb, bin, target_2qb_gate,
                    cx_fidelity, allow_swaps);
              }
              current_interaction[i_vec[i].q0] = -1;
              current_interaction[i_vec[i].q1] = -1;
            }
            if (!is_final_q_type(type)) {
              current_edge_on_qb[q] =
                  circ.get_next_edge(*v, current_edge_on_qb[q]);
            }
          }
        } else if (
            is_gate_type(type) &&
            circ.n_in_edges_of_type(*v, EdgeType::Quantum) == 2) {
          // Check for 2qb gate
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
          } else {
            // End any other interactions on q0
            if (i0 != -1) {
              if (i_vec[i0].count >= 2) {
                // Replace subcircuit
                success |= replace_two_qubit_interaction(
                    circ, i_vec[i0], current_edge_on_qb, bin, target_2qb_gate,
                    cx_fidelity, allow_swaps);
              }
              current_interaction[i_vec[i0].q0] = -1;
              current_interaction[i_vec[i0].q1] = -1;
            }
            // End any other interactions on q1
            if (i1 != -1) {
              if (i_vec[i1].count >= 2) {
                // Replace subcircuit

                success |= replace_two_qubit_interaction(
                    circ, i_vec[i1], current_edge_on_qb, bin, target_2qb_gate,
                    cx_fidelity, allow_swaps);
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
          }
        } else {
          // We don't care about single-qubit vertices, so just update edges
          // and add vertices if interactions exist
          for (port_t i = 0; i < circ.n_in_edges(*v); i++) {
            Qubit q = v_to_qb.at({*v, i});
            current_edge_on_qb[q] =
                circ.get_next_edge(*v, current_edge_on_qb[q]);
            int inter = current_interaction[q];
            if (inter != -1) {
              i_vec[inter].vertices.insert(*v);
            }
          }
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    if (success) {
      squash_1qb_to_tk1().apply(circ);
    }
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
            Subcircuit sub{{in_e}, {out_e}, {}, bin};
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
            sub = Subcircuit{{in_e}, {out_e}, {}, bin};
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

Transform ZZPhase_to_Rz() {
  // basic optimisation, replace ZZPhase with two Rz(1)
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexSet bin;

    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      if (op->get_type() == OpType::ZZPhase) {
        auto params = op->get_params();
        TKET_ASSERT(params.size() == 1);
        // evaluate
        std::optional<double> param = eval_expr(params[0]);
        if (param.has_value() && std::abs(param.value()) == 1.0) {
          success = true;
          // basic optimisation, replace ZZPhase with two Rz(1)
          Circuit replacement(2);
          replacement.add_op<unsigned>(OpType::Rz, 1.0, {0});
          replacement.add_op<unsigned>(OpType::Rz, 1.0, {1});
          circ.substitute(replacement, v, Circuit::VertexDeletion::No);
          bin.insert(v);
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform normalise_TK2() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexSet bin;

    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      bool conditional = false;
      while (op->get_type() == OpType::Conditional) {
        op = static_cast<const Conditional &>(*op).get_op();
        conditional = true;
      }
      if (op->get_type() == OpType::TK2) {
        auto params = op->get_params();
        TKET_ASSERT(params.size() == 3);
        if (!in_weyl_chamber({params[0], params[1], params[2]})) {
          success = true;
          if (conditional) {
            circ.substitute_conditional(
                CircPool::TK2_using_normalised_TK2(
                    params[0], params[1], params[2]),
                v, Circuit::VertexDeletion::No);
          } else {
            circ.substitute(
                CircPool::TK2_using_normalised_TK2(
                    params[0], params[1], params[2]),
                v, Circuit::VertexDeletion::No);
          }
          bin.insert(v);
        }
      }
    }

    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    return success;
  });
}

Transform round_angles(unsigned n, bool only_zeros) {
  if (n >= 32) {
    throw std::invalid_argument("Precision parameter must be less than 32.");
  }
  return Transform([n, only_zeros](Circuit &circ) {
    bool changed = false;
    VertexSet bin;
    const unsigned pow2n = 1 << n;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      bool conditional = false;
      while (op->get_type() == OpType::Conditional) {
        op = static_cast<const Conditional &>(*op).get_op();
        conditional = true;
      }
      if (op->get_desc().is_gate()) {
        std::vector<Expr> params = op->get_params();
        if (!params.empty()) {
          std::vector<Expr> new_params;
          for (const Expr &e : params) {
            std::optional<double> eval = eval_expr(e);
            if (eval) {
              if (only_zeros) {
                new_params.push_back(
                    pow2n * std::abs(eval.value()) < 0.5 ? 0. : e);
              } else {
                new_params.push_back(nearbyint(pow2n * eval.value()) / pow2n);
              }
            } else {
              new_params.push_back(e);
            }
          }
          if (params != new_params) {
            unsigned n_qb = op->n_qubits();
            Circuit replacement(n_qb);
            if (std::any_of(
                    new_params.begin(), new_params.end(),
                    [](const Expr &e) { return (e != 0.); })) {
              std::vector<Qubit> args;
              for (unsigned i = 0; i < n_qb; i++) {
                args.push_back(Qubit(i));
              }
              replacement.add_op(op->get_type(), new_params, args);
            }
            if (conditional) {
              circ.substitute_conditional(
                  replacement, v, Circuit::VertexDeletion::No);
            } else {
              circ.substitute(replacement, v, Circuit::VertexDeletion::No);
            }
            bin.insert(v);
            changed = true;
          }
        }
      }
    }

    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    return changed;
  });
}

}  // namespace tket::Transforms

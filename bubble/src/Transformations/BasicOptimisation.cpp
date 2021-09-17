// Copyright 2019-2021 Cambridge Quantum Computing
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

#include <optional>

#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "Circuit/Command.hpp"
#include "Circuit/DAGDefs.hpp"
#include "Gate/Gate.hpp"
#include "Gate/GatePtr.hpp"
#include "Gate/Rotation.hpp"
#include "Transform.hpp"
#include "Utils/Assert.hpp"
#include "Utils/EigenConfig.hpp"

namespace tket {

static bool redundancy_removal(Circuit &circ);
static bool remove_redundancy(
    Circuit &circ, const Vertex &vert, VertexList &bin,
    std::set<IVertex> &new_affected_verts, IndexMap &im);
static bool commute_singles_to_front(Circuit &circ);
static bool squash_to_pqp(
    Circuit &circ, OpType q, OpType p, bool strict = false);

Transform Transform::remove_redundancies() {
  return Transform(redundancy_removal);
}

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
        z_followed_by_measures &= op->commutes_with_basis(Pauli::Z, port);
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

Transform Transform::squash_1qb_to_tk1() {
  return decompose_ZY() >> squash_1qb_to_pqp(OpType::Ry, OpType::Rz, true) >>
         decompose_ZYZ_to_TK1();
}

Transform Transform::commute_through_multis() {
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
              prev_op->commuting_basis(ports.second);

          if (curr_op->commutes_with_basis(prev_colour, ports.first)) {
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

static bool replace_two_qubit_interaction(
    Circuit &circ, Transform::Interaction &i,
    std::map<Qubit, Edge> &current_edges, VertexList &bin,
    const double cx_fidelity = 1.) {
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
  Circuit replacement = two_qubit_canonical(mat, cx_fidelity);
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

Transform Transform::commute_and_combine_HQS2() {
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
Transform Transform::two_qubit_squash(double cx_fidelity) {
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

Transform Transform::reduce_XZ_chains() {
  return Transform([](Circuit &circ) {
    return squash_to_pqp(circ, OpType::Rx, OpType::Rz);
  });
}

Transform Transform::squash_1qb_to_pqp(
    const OpType &q, const OpType &p, bool strict) {
  return Transform(
      [=](Circuit &circ) { return squash_to_pqp(circ, q, p, strict); });
}

class Squasher {
 public:
  Squasher(Circuit &circ, OpType p, OpType q, bool smart_squash = true)
      : circ_(circ),
        p_(p),
        q_(q),
        success(false),
        bin(),
        smart_squash_(smart_squash) {
    if (!(p == OpType::Rx || p == OpType::Ry || p == OpType::Rz) ||
        !(q == OpType::Rx || q == OpType::Ry || q == OpType::Rz)) {
      throw std::logic_error(
          "Can only reduce chains of single qubit rotations");
    }
    if (p == q) {
      throw std::logic_error(
          "Requires two different bases to perform single qubit "
          "rotations");
    }
  }

  // this squashes the circuit backwards, so that rotations get pushed towards
  // the front, see the design choice notes in confluence
  // https://cqc.atlassian.net/l/c/17xm5hvp
  bool squash() {
    VertexVec outputs = circ_.q_outputs();
    for (VertexVec::const_iterator i = outputs.begin(); i != outputs.end();
         ++i) {
      squash_wire(i);
    }
    circ_.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  }

 private:
  Rotation merge_rotations(
      OpType r, const VertexList &chain,
      VertexList::const_iterator &iter) const {
    Expr total_angle(0);
    while (iter != chain.end()) {
      const Op_ptr rot_op = circ_.get_Op_ptr_from_Vertex(*iter);
      if (rot_op->get_type() != r) {
        break;
      }
      total_angle += rot_op->get_params()[0];
      iter++;
    }
    return Rotation(r, total_angle);
  }

  static bool fixup_angles(Expr &angle_p1, Expr &angle_q, Expr &angle_p2) {
    if (equiv_val(angle_q, 1., 2) && !equiv_0(angle_p2, 4)) {
      // Prefer --P(p1-p2)--Q(...)--P(0)--
      // Only occurs if angle_q is pi or 3pi and angle_p2 is non-zero
      angle_p1 -= angle_p2;
      angle_p2 = 0;
      return true;
    } else if (equiv_val(angle_p2, 1., 4)) {
      // Then prefer --P(p1+p2)--Q(-q)--P(0)--
      // Only occurs if angle_p2 is pi
      angle_p1 += 1;
      angle_q *= -1;
      angle_p2 = 0;
      return true;
    } else if (equiv_val(angle_p2, 3., 4)) {
      // Then prefer --P(p1+p2)--Q(-q)--P(0)--
      // Only occurs if angle_p2 is 3pi
      angle_p1 += 3;
      angle_q *= -1;
      angle_p2 = 0;
      return true;
    } else if (equiv_val(angle_p1, 1., 4) && !equiv_0(angle_p2, 4)) {
      // Then prefer --P(0)--Q(-q)--P(p1+p2)--
      // Only occurs if angle_p1 is pi and angle_p2 is non-zero
      angle_q *= -1;
      angle_p2 += 1;
      angle_p1 = 0;
      return true;
    } else if (equiv_val(angle_p1, 3., 4) && !equiv_0(angle_p2, 4)) {
      // Then prefer --P(0)--Q(-q)--P(p1+p2)--
      // Only occurs if angle_p1 is 3pi and angle_p2 is non-zero
      angle_q *= -1;
      angle_p2 += 3;
      angle_p1 = 0;
      return true;
    } else
      return false;
  }

  bool is_same_chain(
      const Circuit &replacement, const VertexList &rotation_chain) const {
    auto orig_rot = rotation_chain.begin();
    if (rotation_chain.size() != replacement.n_gates()) {
      return false;
    }
    BGL_FORALL_VERTICES(new_rot, replacement.dag, DAG) {
      Op_ptr new_rot_op = replacement.get_Op_ptr_from_Vertex(new_rot);
      if (!is_boundary_q_type(new_rot_op->get_type())) {
        Op_ptr orig_rot_op = circ_.get_Op_ptr_from_Vertex(*orig_rot);
        if (!(*new_rot_op == *orig_rot_op)) {
          return false;
        }
        ++orig_rot;
      }
    }
    return true;
  }

  std::tuple<Expr, Expr, Expr> pqp_from_chain(
      const VertexList &rotation_chain, bool invert_pqp = false) const {
    OpType p = (invert_pqp) ? q_ : p_;
    OpType q = (invert_pqp) ? p_ : q_;

    // Construct list of merged rotations
    std::list<Rotation> rots;
    VertexList::const_iterator iter = rotation_chain.begin();
    while (iter != rotation_chain.end()) {
      // Merge next q rotations
      rots.push_back(merge_rotations(q, rotation_chain, iter));
      // Merge next p rotations
      rots.push_back(merge_rotations(p, rotation_chain, iter));
    }

    // Perform any cancellations
    std::list<Rotation>::iterator r = rots.begin();
    while (r != rots.end()) {
      if (r->is_id()) {
        r = rots.erase(r);
        if (r != rots.begin() && r != rots.end()) {
          std::prev(r)->apply(*r);
          r = rots.erase(r);
          r--;
        }
      } else
        r++;
    }

    // Extract any P rotations from the beginning and end of the list
    Expr p1 = 0, p2 = 0;
    std::list<Rotation>::iterator i1 = rots.begin();
    if (i1 != rots.end()) {
      std::optional<Expr> a = i1->angle(p);
      if (a) {
        p1 = a.value();
        rots.pop_front();
      }
    }
    std::list<Rotation>::reverse_iterator i2 = rots.rbegin();
    if (i2 != rots.rend()) {
      std::optional<Expr> a = i2->angle(p);
      if (a) {
        p2 = a.value();
        rots.pop_back();
      }
    }

    // Finish up:
    Rotation R = {};
    for (auto rot : rots) {
      R.apply(rot);
    }
    std::tuple<Expr, Expr, Expr> pqp = R.to_pqp(p, q);
    std::get<0>(pqp) += p1;
    std::get<2>(pqp) += p2;
    return pqp;
  }

  void squash_rotations(const VertexList &rotation_chain) {
    // TODO:: break chain up with classical control

    // smart squashing: choose pqp (default) or qpq depending on next gate
    bool choose_qpq = false;
    // smart squashing: flag if last rotation can be commuted through next gate
    bool commute_through = false;

    const Op_ptr next_op = circ_.get_Op_ptr_from_Vertex(v);
    const OpType next_op_type = next_op->get_type();
    if (smart_squash_ && is_gate_type(next_op_type)) {
      const port_t source_port = circ_.get_source_port(e);
      const std::optional<Pauli> commutation_colour =
          next_op->commuting_basis(source_port);

      Gate P(p_, {0}, 1);
      Gate Q(q_, {0}, 1);
      if (P.commutes_with_basis(commutation_colour, 0)) {
        commute_through = true;
      } else if (Q.commutes_with_basis(commutation_colour, 0)) {
        choose_qpq = true;
        commute_through = true;
      }
    }

    // swap p, q if required
    OpType p = choose_qpq ? q_ : p_;
    OpType q = choose_qpq ? p_ : q_;

    auto angles = pqp_from_chain(rotation_chain, choose_qpq);

    Expr angle_p1 = std::get<0>(angles);
    Expr angle_q = std::get<1>(angles);
    Expr angle_p2 = std::get<2>(angles);
    (void)fixup_angles(angle_p1, angle_q, angle_p2);

    Circuit replacement(1);
    if (!commute_through) {
      replacement.add_op<unsigned>(p, angle_p1, {0});
    }
    replacement.add_op<unsigned>(q, angle_q, {0});
    replacement.add_op<unsigned>(p, angle_p2, {0});
    redundancy_removal(replacement);

    // check if replacement is any different from original chain
    if (!is_same_chain(replacement, rotation_chain)) {
      success = true;

      // replace with new rotations in circuit
      Subcircuit sub = {
          {e}, {circ_.get_nth_out_edge(rotation_chain.back(), 0)}};
      port_t port = circ_.get_source_port(e);
      circ_.substitute(replacement, sub, Circuit::VertexDeletion::No);
      e = circ_.get_nth_out_edge(v, port);
      bin.insert(bin.end(), rotation_chain.begin(), rotation_chain.end());
      // add gate commuted through v
      if (commute_through) {
        Edge last_e = circ_.get_last_edge(v, e);
        Subcircuit before_v{{last_e}, {last_e}};
        Circuit leftover_p_gate(1);
        leftover_p_gate.add_op<unsigned>(p, angle_p1, {0});
        circ_.substitute(
            leftover_p_gate, before_v, Circuit::VertexDeletion::No);
      }
    }
  }

  void squash_wire(VertexVec::const_iterator i) {
    e = circ_.get_nth_in_edge(*i, 0);
    v = circ_.source(e);
    VertexList rotation_chain;
    while (true) {
      OpType v_type = circ_.get_OpType_from_Vertex(v);
      if (v_type == p_ || v_type == q_) {
        rotation_chain.push_front(v);
      } else if (!rotation_chain.empty()) {
        squash_rotations(rotation_chain);
        rotation_chain.clear();
      }
      if (is_initial_q_type(v_type)) {
        break;
      }
      e = circ_.get_last_edge(v, e);
      v = circ_.source(e);
    }
  }

  Circuit &circ_;
  OpType p_;
  OpType q_;
  bool success;
  VertexList bin;
  bool smart_squash_;
  Edge e;
  Vertex v;
};

static bool squash_to_pqp(Circuit &circ, OpType q, OpType p, bool strict) {
  Squasher s(circ, p, q, !strict);
  return s.squash();
}

static bool standard_squash(
    Circuit &circ, const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement) {
  bool success = false;
  for (OpType ot : singleqs) {
    if (!is_single_qubit_type(ot))
      throw NotValid(
          "OpType given to standard_squash is not a single qubit "
          "gate");
  }
  for (const Vertex in : circ.q_inputs()) {
    Vertex v = in;
    port_t port = 0;
    VertexList single_chain;
    Rotation combined;
    std::optional<std::pair<std::list<VertPort>, unsigned>> condition;
    while (true) {
      Op_ptr v_op = circ.get_Op_ptr_from_Vertex(v);
      OpType v_type = v_op->get_type();
      bool squash_chain = false;
      bool add_to_chain = false;
      std::optional<std::pair<std::list<VertPort>, unsigned>> this_condition =
          std::nullopt;
      if (v_type == OpType::Conditional) {
        const Conditional &cond_op = static_cast<const Conditional &>(*v_op);
        EdgeVec ins = circ.get_in_edges(v);
        this_condition = std::pair<std::list<VertPort>, unsigned>();
        for (port_t p = 0; p < cond_op.get_width(); ++p) {
          Edge in_p = ins.at(p);
          VertPort vp = {circ.source(in_p), circ.get_source_port(in_p)};
          this_condition->first.push_back(vp);
        }
        this_condition->second = cond_op.get_value();
        v_op = cond_op.get_op();
        v_type = v_op->get_type();
      }
      if (condition != this_condition) squash_chain = true;
      unsigned q_ins = circ.n_in_edges_of_type(v, EdgeType::Quantum);
      if ((q_ins != 1) | (singleqs.find(v_type) == singleqs.end()) |
          is_projective_type(v_type))
        squash_chain = true;
      else
        add_to_chain = true;
      if (squash_chain && !single_chain.empty()) {
        auto [a, b, c] = combined.to_pqp(OpType::Rz, OpType::Rx);
        Circuit replacement = tk1_replacement(c, b, a);
        if (replacement.n_gates() < single_chain.size()) {
          BGL_FORALL_VERTICES(rv, replacement.dag, DAG) {
            OpType v_type = replacement.get_OpType_from_Vertex(rv);
            if (!is_boundary_q_type(v_type) &&
                singleqs.find(v_type) == singleqs.end())
              throw NotValid(
                  "tk1_replacement given to standard_squash "
                  "does not preserve gate set");
          }
          if (condition) {
            circ.substitute_conditional(
                replacement, single_chain.front(), Circuit::VertexDeletion::No);
          } else {
            circ.substitute(
                replacement, single_chain.front(), Circuit::VertexDeletion::No);
          }
          circ.remove_vertices(
              single_chain, Circuit::GraphRewiring::Yes,
              Circuit::VertexDeletion::Yes);
          success = true;
        }
        single_chain = VertexList();
        combined = Rotation();
        condition = std::nullopt;
      }
      if (is_final_q_type(v_type)) break;
      if (add_to_chain) {
        if (single_chain.empty()) condition = this_condition;
        single_chain.push_back(v);
        std::vector<Expr> angs = as_gate_ptr(v_op)->get_tk1_angles();
        combined.apply(Rotation(OpType::Rz, angs.at(2)));
        combined.apply(Rotation(OpType::Rx, angs.at(1)));
        combined.apply(Rotation(OpType::Rz, angs.at(0)));
      }
      Edge e = circ.get_nth_out_edge(v, port);
      v = circ.target(e);
      port = circ.get_target_port(e);
    }
  }
  return success;
}

Transform Transform::squash_factory(
    const OpTypeSet &singleqs,
    const std::function<Circuit(const Expr &, const Expr &, const Expr &)>
        &tk1_replacement) {
  return Transform([=](Circuit &circ) {
    return standard_squash(circ, singleqs, tk1_replacement);
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
Transform Transform::commute_SQ_gates_through_SWAPS(
    const avg_node_errors_t &node_errors) {
  return commute_SQ_gates_through_SWAPS_helper(
      DeviceCharacterisation(node_errors));
}
Transform Transform::commute_SQ_gates_through_SWAPS(
    const op_node_errors_t &node_errors) {
  return commute_SQ_gates_through_SWAPS_helper(
      DeviceCharacterisation(node_errors));
}

}  // namespace tket

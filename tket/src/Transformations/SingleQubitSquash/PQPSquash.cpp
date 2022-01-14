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

#include "Gate/Rotation.hpp"
#include "SingleQubitSquash.hpp"

namespace tket {

static bool fixup_angles(
    Expr &angle_p1, Expr &angle_q, Expr &angle_p2, bool reversed = false);
static bool redundancy_removal(Circuit &circ);

/**
 * @brief Implements the Squasher interface for SingleQubitSquash
 *
 * The PQP Squasher squashes chains of single qubit gates to minimal sequences
 * of any two type of rotations (Rx, Ry or Rz), using Euler angle
 * decompositions.
 *
 * Decompositions can either be strict or smart. Strict decompositions
 * will always decompose to P-Q-P, whereas smart will decompose to P-Q-P
 * or Q-P-Q and will always try to commute P or Q through multi-qubit gates.
 * For smart decomposition, swapping P and Q gives the same result.
 *
 * Note that even in strict mode, if a rotation in the PQP decomposition
 * has angle 0, it will be omitted.
 *
 * The squash is made backwards, so that rotations get pushed towards
 * the front, see the design choice notes in confluence
 * https://cqc.atlassian.net/l/c/17xm5hvp
 *
 * The PQPSquasher shouldn't have to know if the squash is made forwards
 * or backwards. Here, it is used only in fixup_angles, because the convention
 * is different (for consistency to the user)
 */
class PQPSquasher {
 public:
  PQPSquasher(
      OpType p, OpType q, bool smart_squash = true, bool reversed = false)
      : p_(p),
        q_(q),
        smart_squash_(smart_squash),
        reversed_(reversed),
        rotation_chain() {
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

  void clear() { rotation_chain.clear(); }

  bool accepts(OpType type) const { return type == p_ || type == q_; }

  void append(Gate_ptr gp) {
    if (!accepts(gp->get_type())) {
      throw NotValid("PQPSquasher: cannot append OpType");
    }
    rotation_chain.push_back(gp);
  }

  std::pair<Circuit, Gate_ptr> flush(
      std::optional<Pauli> commutation_colour = std::nullopt) const {
    bool commute_through = false;
    OpType p = p_, q = q_;

    if (smart_squash_ && commutation_colour.has_value()) {
      Gate P(p_, {0}, 1);
      Gate Q(q_, {0}, 1);
      if (P.commutes_with_basis(commutation_colour, 0)) {
        commute_through = true;
      } else if (Q.commutes_with_basis(commutation_colour, 0)) {
        commute_through = true;
        p = q_;
        q = p_;
      }
    }

    // Construct list of merged rotations
    std::list<Rotation> rots;
    auto iter = rotation_chain.cbegin();
    while (iter != rotation_chain.cend()) {
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
    std::tuple<Expr, Expr, Expr> angles = R.to_pqp(p, q);

    Expr angle_p1 = std::get<0>(angles) + p1;
    Expr angle_q = std::get<1>(angles);
    Expr angle_p2 = std::get<2>(angles) + p2;
    fixup_angles(angle_p1, angle_q, angle_p2, reversed_);

    Circuit replacement(1);
    Gate_ptr left_over_gate = nullptr;
    if (!equiv_0(angle_p1, 4)) {
      if (equiv_0(angle_q, 4) && equiv_0(angle_p2, 4) && commute_through) {
        left_over_gate =
            std::make_shared<Gate>(p, std::vector<Expr>{angle_p1}, 1);
      } else {
        replacement.add_op<unsigned>(p, angle_p1, {0});
      }
    }
    if (!equiv_0(angle_q, 4)) {
      replacement.add_op<unsigned>(q, angle_q, {0});
    }
    if (!equiv_0(angle_p2, 4)) {
      if (commute_through) {
        left_over_gate =
            std::make_shared<Gate>(p, std::vector<Expr>{angle_p2}, 1);
      } else {
        replacement.add_op<unsigned>(p, angle_p2, {0});
      }
    }
    redundancy_removal(replacement);
    return {replacement, left_over_gate};
  }

 private:
  static Rotation merge_rotations(
      OpType r, const std::vector<Gate_ptr> &chain,
      std::vector<Gate_ptr>::const_iterator &iter) {
    Expr total_angle(0);
    while (iter != chain.end()) {
      const Gate_ptr rot_op = *iter;
      if (rot_op->get_type() != r) {
        break;
      }
      total_angle += rot_op->get_params()[0];
      iter++;
    }
    return Rotation(r, total_angle);
  }

  OpType p_;
  OpType q_;
  bool smart_squash_;
  bool reversed_;
  std::vector<Gate_ptr> rotation_chain;
};

static bool squash_to_pqp(
    Circuit &circ, OpType q, OpType p, bool strict = false) {
  bool reverse = true;
  PQPSquasher squasher(p, q, !strict, reverse);
  return SingleQubitSquash(squasher, reverse).squash(circ);
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

static bool fixup_angles(
    Expr &angle_p1, Expr &angle_q, Expr &angle_p2, bool reversed) {
  bool success = false;
  if (reversed) {
    std::swap(angle_p1, angle_p2);
    angle_p1 *= -1;
    angle_q *= -1;
    angle_p2 *= -1;
  }
  if (equiv_val(angle_q, 1., 2) && !equiv_0(angle_p2, 4)) {
    // Prefer --P(p1-p2)--Q(...)--P(0)--
    // Only occurs if angle_q is pi or 3pi and angle_p2 is non-zero
    angle_p1 -= angle_p2;
    angle_p2 = 0;
    success = true;
  } else if (equiv_val(angle_p2, 1., 4)) {
    // Then prefer --P(p1+p2)--Q(-q)--P(0)--
    // Only occurs if angle_p2 is pi
    angle_p1 += 1;
    angle_q *= -1;
    angle_p2 = 0;
    success = true;
  } else if (equiv_val(angle_p2, 3., 4)) {
    // Then prefer --P(p1+p2)--Q(-q)--P(0)--
    // Only occurs if angle_p2 is 3pi
    angle_p1 += 3;
    angle_q *= -1;
    angle_p2 = 0;
    success = true;
  } else if (equiv_val(angle_p1, 1., 4) && !equiv_0(angle_p2, 4)) {
    // Then prefer --P(0)--Q(-q)--P(p1+p2)--
    // Only occurs if angle_p1 is pi and angle_p2 is non-zero
    angle_q *= -1;
    angle_p2 += 1;
    angle_p1 = 0;
    success = true;
  } else if (equiv_val(angle_p1, 3., 4) && !equiv_0(angle_p2, 4)) {
    // Then prefer --P(0)--Q(-q)--P(p1+p2)--
    // Only occurs if angle_p1 is 3pi and angle_p2 is non-zero
    angle_q *= -1;
    angle_p2 += 3;
    angle_p1 = 0;
    success = true;
  }
  if (reversed) {
    std::swap(angle_p1, angle_p2);
    angle_p1 *= -1;
    angle_q *= -1;
    angle_p2 *= -1;
  }
  return success;
}

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

}  // namespace tket
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

#include "tket/Transformations/PQPSquash.hpp"

#include <memory>

#include "tket/Circuit/DAGDefs.hpp"
#include "tket/Gate/Rotation.hpp"
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/Transform.hpp"
#include "tket/Utils/Expression.hpp"

namespace tket::Transforms {

static bool fixup_angles(
    Expr &angle_p1, Expr &angle_q, Expr &angle_p2, bool reversed = false);
static Rotation merge_rotations(
    OpType r, const std::vector<Gate_ptr> &chain,
    std::vector<Gate_ptr>::const_iterator &iter);

PQPSquasher::PQPSquasher(OpType p, OpType q, bool smart_squash, bool reversed)
    : p_(p),
      q_(q),
      smart_squash_(smart_squash),
      reversed_(reversed),
      rotation_chain() {
  if (!(p == OpType::Rx || p == OpType::Ry || p == OpType::Rz) ||
      !(q == OpType::Rx || q == OpType::Ry || q == OpType::Rz)) {
    throw std::logic_error("Can only reduce chains of single qubit rotations");
  }
  if (p == q) {
    throw std::logic_error(
        "Requires two different bases to perform single qubit "
        "rotations");
  }
}

bool PQPSquasher::accepts(Gate_ptr gp) const {
  OpType type = gp->get_type();
  return type == p_ || type == q_;
}

void PQPSquasher::append(Gate_ptr gp) {
  if (!accepts(gp)) {
    throw BadOpType("PQPSquasher: cannot append OpType", gp->get_type());
  }
  rotation_chain.push_back(gp);
}

std::pair<Circuit, Gate_ptr> PQPSquasher::flush(
    std::optional<Pauli> commutation_colour) const {
  bool commute_through = false;
  OpType p = p_, q = q_;

  if (smart_squash_ && commutation_colour.has_value()) {
    // Using an arbitrary non-zero angle to obtain the commutation for p_/q_.
    Gate P(p_, {0.123}, 1);
    Gate Q(q_, {0.123}, 1);
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

void PQPSquasher::clear() { rotation_chain.clear(); }

std::unique_ptr<AbstractSquasher> PQPSquasher::clone() const {
  return std::make_unique<PQPSquasher>(*this);
}

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

static bool squash_to_pqp(
    Circuit &circ, OpType q, OpType p, bool strict = false) {
  bool reverse = true;
  auto squasher = std::make_unique<PQPSquasher>(p, q, !strict, reverse);
  return SingleQubitSquash(std::move(squasher), circ, reverse).squash();
}

Transform reduce_XZ_chains() {
  return Transform([](Circuit &circ) {
    return squash_to_pqp(circ, OpType::Rx, OpType::Rz);
  });
}

Transform squash_1qb_to_pqp(const OpType &q, const OpType &p, bool strict) {
  return Transform(
      [=](Circuit &circ) { return squash_to_pqp(circ, q, p, strict); });
}

// To squash to TK1:
// - we first decompose to ZYZ. Doing this was found to reduce the size of
//   symbolic expressions
// - we then redecompose to ZXZ, so that we can commute Rz or Rx rotation past
//   multi-qubit gates (most usual multi-qb gates commute with X or Z)
// - Rz and Rx rotations can then be straight-forwardly combined into TK1s.
Transform squash_1qb_to_tk1() {
  return Transforms::decompose_ZY() >>
         squash_1qb_to_pqp(OpType::Ry, OpType::Rz, true) >>
         Transforms::decompose_ZX() >>
         squash_1qb_to_pqp(OpType::Rx, OpType::Rz, true) >>
         Transforms::decompose_ZXZ_to_TK1();
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

}  // namespace tket::Transforms

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

#include "Circuit/CircPool.hpp"
#include "Circuit/CircUtils.hpp"
#include "Gate/GatePtr.hpp"
#include "Replacement.hpp"
#include "Transform.hpp"
#include "Utils/TketLog.hpp"

namespace tket {

static bool standard_rebase(
    Circuit& circ, const OpTypeSet& multiqs, const Circuit& cx_replacement,
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement);

Transform Transform::rebase_factory(
    const OpTypeSet& multiqs, const Circuit& cx_replacement,
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  return Transform([=](Circuit& circ) {
    return standard_rebase(
        circ, multiqs, cx_replacement, singleqs, tk1_replacement);
  });
}

static bool standard_rebase(
    Circuit& circ, const OpTypeSet& multiqs, const Circuit& cx_replacement,
    const OpTypeSet& singleqs,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  bool success = false;
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    unsigned n_qubits = circ.n_in_edges_of_type(v, EdgeType::Quantum);
    if (n_qubits <= 1) continue;
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (multiqs.find(type) != multiqs.end() || type == OpType::CX ||
        type == OpType::Barrier)
      continue;
    // need to convert
    Circuit replacement = CX_circ_from_multiq(op);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    bin.push_back(v);
    success = true;
  }
  if (multiqs.find(OpType::CX) == multiqs.end()) {
    const Op_ptr cx_op = get_op_ptr(OpType::CX);
    success = circ.substitute_all(cx_replacement, cx_op) | success;
  }
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges_of_type(v, EdgeType::Quantum) != 1 ||
        circ.n_in_edges_of_type(v, EdgeType::Quantum) != 1)
      continue;
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (!is_gate_type(type) || is_projective_type(type) ||
        singleqs.find(type) != singleqs.end())
      continue;
    // need to convert
    std::vector<Expr> tk1_angles = as_gate_ptr(op)->get_tk1_angles();
    Circuit replacement =
        tk1_replacement(tk1_angles[0], tk1_angles[1], tk1_angles[2]);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    circ.add_phase(tk1_angles[3]);
    bin.push_back(v);
    success = true;
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

Transform Transform::rebase_tket() {
  std::function<Circuit(const Expr&, const Expr&, const Expr&)> tk1_to_tk1 =
      [](const Expr& alpha, const Expr& beta, const Expr& gamma) {
        Circuit c(1);
        c.add_op<unsigned>(OpType::TK1, {alpha, beta, gamma}, {0});
        return c;
      };
  return rebase_factory(
      {OpType::CX}, CircPool::CX(), {OpType::TK1}, tk1_to_tk1);
}

Circuit Transform::tk1_to_PhasedXRz(
    const Expr& alpha, const Expr& beta, const Expr& gamma) {
  Circuit c(1);
  if (equiv_expr(beta, 1)) {
    // Angles β ∈ {π, 3π}
    c.add_op<unsigned>(OpType::PhasedX, {beta, (alpha - gamma) / 2.}, {0});
  } else if (equiv_expr(beta, 0)) {
    // Angle β ∈ {0, 2π}
    c.add_op<unsigned>(OpType::Rz, alpha + beta + gamma, {0});
  } else {
    c.add_op<unsigned>(OpType::Rz, alpha + gamma, {0});
    c.add_op<unsigned>(OpType::PhasedX, {beta, alpha}, {0});
  }
  Transform::remove_redundancies().apply(c);
  return c;
}

Transform Transform::rebase_cirq() {
  return rebase_factory(
      {OpType::CZ}, CircPool::H_CZ_H(), {OpType::PhasedX, OpType::Rz},
      Transform::tk1_to_PhasedXRz);
}

Circuit Transform::tk1_to_rzrx(
    const Expr& alpha, const Expr& beta, const Expr& gamma) {
  Circuit c(1);
  c.add_op<unsigned>(OpType::Rz, gamma, {0});
  c.add_op<unsigned>(OpType::Rx, beta, {0});
  c.add_op<unsigned>(OpType::Rz, alpha, {0});
  Transform::remove_redundancies().apply(c);
  return c;
}

Circuit Transform::tk1_to_rzh(
    const Expr& alpha, const Expr& beta, const Expr& gamma) {
  Circuit c(1);
  std::optional<unsigned> cliff = equiv_Clifford(beta, 4);
  if (cliff) {
    switch (*cliff % 4) {
      case 0: {
        c.add_op<unsigned>(OpType::Rz, gamma + alpha, {0});
        break;
      }
      case 1: {
        c.add_op<unsigned>(OpType::Rz, gamma - 0.5, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, alpha - 0.5, {0});
        c.add_phase(-0.5);
        break;
      }
      case 2: {
        c.add_op<unsigned>(OpType::Rz, gamma - alpha, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, 1., {0});
        c.add_op<unsigned>(OpType::H, {0});
        break;
      }
      case 3: {
        c.add_op<unsigned>(OpType::Rz, gamma + 0.5, {0});
        c.add_op<unsigned>(OpType::H, {0});
        c.add_op<unsigned>(OpType::Rz, alpha + 0.5, {0});
        c.add_phase(-0.5);
        break;
      }
    }
    if (cliff >= 4u) c.add_phase(1.);
  } else {
    c.add_op<unsigned>(OpType::Rz, gamma, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, beta, {0});
    c.add_op<unsigned>(OpType::H, {0});
    c.add_op<unsigned>(OpType::Rz, alpha, {0});
  }
  Transform::remove_redundancies().apply(c);
  return c;
}

static unsigned int_half(const Expr& angle) {
  // Assume angle is an even integer
  double eval = eval_expr(angle).value();
  return lround(eval / 2);
}

Circuit Transform::tk1_to_rzsx(
    const Expr& alpha, const Expr& beta, const Expr& gamma) {
  Circuit c(1);
  Expr correction_phase = 0;
  if (equiv_0(beta)) {
    // b = 2k, if k is odd, then Rx(b) = -I
    c.add_op<unsigned>(OpType::Rz, alpha + gamma, {0});
    correction_phase = int_half(beta);
  } else if (equiv_0(beta + 1)) {
    // Use Rx(2k-1) = i(-1)^{k}SxSx
    correction_phase = -0.5 + int_half(beta - 1);
    if (equiv_0(alpha - gamma)) {
      // a - c = 2m
      // overall operation is (-1)^{m}Rx(2k -1)
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      correction_phase += int_half(alpha - gamma);
    } else {
      c.add_op<unsigned>(OpType::Rz, gamma, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::SX, {0});
      c.add_op<unsigned>(OpType::Rz, alpha, {0});
    }
  } else if (equiv_0(beta - 0.5) && equiv_0(alpha) && equiv_0(gamma)) {
    // a = 2k, b = 2m+0.5, c = 2n
    // Rz(2k)Rx(2m + 0.5)Rz(2n) = (-1)^{k+m+n}e^{-i \pi /4} SX
    c.add_op<unsigned>(OpType::SX, {0});
    correction_phase =
        int_half(beta - 0.5) + int_half(alpha) + int_half(gamma) - 0.25;
  } else if (equiv_0(alpha - 0.5) && equiv_0(gamma - 0.5)) {
    // Rz(2k + 0.5)Rx(b)Rz(2m + 0.5) = -i(-1)^{k+m}SX.Rz(1-b).SX
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, 1 - beta, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    correction_phase = int_half(alpha - 0.5) + int_half(gamma - 0.5) - 0.5;
  } else {
    c.add_op<unsigned>(OpType::Rz, gamma + 0.5, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, beta - 1, {0});
    c.add_op<unsigned>(OpType::SX, {0});
    c.add_op<unsigned>(OpType::Rz, alpha + 0.5, {0});
    correction_phase = -0.5;
  }
  c.add_phase(correction_phase);
  Transform::remove_redundancies().apply(c);
  return c;
}

Circuit Transform::tk1_to_tk1(
    const Expr& alpha, const Expr& beta, const Expr& gamma) {
  Circuit c(1);
  c.add_op<unsigned>(OpType::TK1, {alpha, beta, gamma}, {0});
  return c;
}

Transform Transform::rebase_HQS() {
  return rebase_factory(
      {OpType::ZZMax}, CircPool::CX_using_ZZMax(),
      {OpType::PhasedX, OpType::Rz}, Transform::tk1_to_PhasedXRz);
}

Transform Transform::rebase_UMD() {
  return rebase_factory(
      {OpType::XXPhase}, CircPool::CX_using_XXPhase_0(),
      {OpType::PhasedX, OpType::Rz}, Transform::tk1_to_PhasedXRz);
}

Transform Transform::rebase_quil() {
  return rebase_factory(
      {OpType::CZ}, CircPool::H_CZ_H(), {OpType::Rx, OpType::Rz}, tk1_to_rzrx);
}

Transform Transform::rebase_pyzx() {
  OpTypeSet pyzx_multiqs = {OpType::SWAP, OpType::CX, OpType::CZ};
  OpTypeSet pyzx_singleqs = {OpType::H, OpType::X,  OpType::Z, OpType::S,
                             OpType::T, OpType::Rx, OpType::Rz};
  return rebase_factory(
      pyzx_multiqs, CircPool::CX(), pyzx_singleqs, tk1_to_rzrx);
}

Transform Transform::rebase_projectq() {
  OpTypeSet projectq_multiqs = {
      OpType::SWAP, OpType::CRz, OpType::CX, OpType::CZ};
  OpTypeSet projectq_singleqs = {OpType::H,  OpType::X, OpType::Y, OpType::Z,
                                 OpType::S,  OpType::T, OpType::V, OpType::Rx,
                                 OpType::Ry, OpType::Rz};
  return rebase_factory(
      projectq_multiqs, CircPool::CX(), projectq_singleqs, tk1_to_rzrx);
}

Transform Transform::rebase_UFR() {
  return rebase_factory(
      {OpType::CX}, CircPool::CX(), {OpType::Rz, OpType::H}, tk1_to_rzh);
}

Transform Transform::rebase_OQC() {
  return rebase_factory(
      {OpType::ECR}, CircPool::CX_using_ECR(), {OpType::Rz, OpType::SX},
      Transform::tk1_to_rzsx);
}

}  // namespace tket

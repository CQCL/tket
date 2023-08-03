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

#include "tket/Transformations/Rebase.hpp"

#include <tkassert/Assert.hpp>
#include <tklog/TketLog.hpp>

#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Gate/GatePtr.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/OpType/OpTypeFunctions.hpp"
#include "tket/Ops/Op.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Replacement.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {

namespace Transforms {

// Rebase a 0- or 1-qubit unitary operation
static Circuit rebase_op(
    Gate_ptr op,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  OpType type = op->get_type();
  if (type == OpType::Phase) {
    Circuit replacement(0);
    replacement.add_phase(op->get_params()[0]);
    return replacement;
  } else {
    std::vector<Expr> tk1_angles = op->get_tk1_angles();
    Circuit replacement =
        tk1_replacement(tk1_angles[0], tk1_angles[1], tk1_angles[2]);
    remove_redundancies().apply(replacement);
    replacement.add_phase(tk1_angles[3]);
    return replacement;
  }
}

static bool standard_rebase(
    Circuit& circ, const OpTypeSet& allowed_gates,
    const Circuit& cx_replacement,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  bool success = false;
  VertexSet bin;
  for (const Vertex& v : circ.all_vertices()) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    unsigned n_qubits = circ.n_in_edges_of_type(v, EdgeType::Quantum);
    if (n_qubits <= 1) continue;
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (allowed_gates.find(type) != allowed_gates.end() || type == OpType::CX ||
        type == OpType::Barrier)
      continue;
    // need to convert
    Circuit replacement = CX_circ_from_multiq(op);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    bin.insert(v);
    success = true;
  }
  if (allowed_gates.find(OpType::CX) == allowed_gates.end()) {
    const Op_ptr cx_op = get_op_ptr(OpType::CX);
    success = circ.substitute_all(cx_replacement, cx_op) | success;
  }
  for (const Vertex& v : circ.all_vertices()) {
    if (bin.contains(v)) continue;
    if (circ.n_in_edges_of_type(v, EdgeType::Quantum) > 1) continue;
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (!is_gate_type(type) || is_projective_type(type) ||
        allowed_gates.find(type) != allowed_gates.end())
      continue;
    // need to convert
    Circuit replacement = rebase_op(as_gate_ptr(op), tk1_replacement);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    bin.insert(v);
    success = true;
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

static bool standard_rebase_via_tk2(
    Circuit& circ, const OpTypeSet& allowed_gates,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk2_replacement) {
  bool success = false;
  VertexSet bin;
  // 1. Replace all multi-qubit gates outside the target gateset to TK2.
  for (const Vertex& v : circ.all_vertices()) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    unsigned n_qubits = circ.n_in_edges_of_type(v, EdgeType::Quantum);
    if (n_qubits <= 1) continue;
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (allowed_gates.contains(type) || type == OpType::TK2 ||
        type == OpType::Barrier)
      continue;
    // need to convert
    Circuit replacement = TK2_circ_from_multiq(op);
    // Find replacement Circuit for all TK2 gates
    VertexSet TK2_bin;
    for (const Vertex& u : replacement.all_vertices()) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(u);
      TKET_ASSERT(op->get_type() != OpType::Conditional);
      if (op->get_type() == OpType::TK2) {
        std::vector<Expr> params = op->get_params();
        TKET_ASSERT(params.size() == 3);
        Circuit u_replacement =
            tk2_replacement(params[0], params[1], params[2]);
        replacement.substitute(u_replacement, u, Circuit::VertexDeletion::No);
        TK2_bin.insert(u);
      }
    }
    Transforms::squash_1qb_to_tk1().apply(replacement);
    remove_redundancies().apply(replacement);
    replacement.remove_vertices(
        TK2_bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    bin.insert(v);
    success = true;
  }

  // 2. Replace 0- and 1-qubit gates by converting to TK1 and replacing.
  for (const Vertex& v : circ.all_vertices()) {
    if (bin.contains(v)) continue;
    if (circ.n_in_edges_of_type(v, EdgeType::Quantum) > 1) continue;
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    bool conditional = op->get_type() == OpType::Conditional;
    if (conditional) {
      const Conditional& cond = static_cast<const Conditional&>(*op);
      op = cond.get_op();
    }
    OpType type = op->get_type();
    if (!is_gate_type(type) || is_projective_type(type) ||
        allowed_gates.contains(type))
      continue;
    // need to convert
    Circuit replacement = rebase_op(as_gate_ptr(op), tk1_replacement);
    if (conditional) {
      circ.substitute_conditional(replacement, v, Circuit::VertexDeletion::No);
    } else {
      circ.substitute(replacement, v, Circuit::VertexDeletion::No);
    }
    bin.insert(v);
    success = true;
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

Transform rebase_factory(
    const OpTypeSet& allowed_gates, const Circuit& cx_replacement,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement) {
  return Transform([=](Circuit& circ) {
    return standard_rebase(
        circ, allowed_gates, cx_replacement, tk1_replacement);
  });
}

Transform rebase_factory_via_tk2(
    const OpTypeSet& allowed_gates,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk1_replacement,
    const std::function<Circuit(const Expr&, const Expr&, const Expr&)>&
        tk2_replacement) {
  return Transform([=](Circuit& circ) {
    return standard_rebase_via_tk2(
        circ, allowed_gates, tk1_replacement, tk2_replacement);
  });
}

Transform rebase_tket() {
  std::function<Circuit(const Expr&, const Expr&, const Expr&)> tk1_to_tk1 =
      [](const Expr& alpha, const Expr& beta, const Expr& gamma) {
        Circuit c(1);
        c.add_op<unsigned>(OpType::TK1, {alpha, beta, gamma}, {0});
        return c;
      };
  return rebase_factory({OpType::CX, OpType::TK1}, CircPool::CX(), tk1_to_tk1);
}

Transform rebase_cirq() {
  return rebase_factory(
      {OpType::CZ, OpType::PhasedX, OpType::Rz}, CircPool::H_CZ_H(),
      CircPool::tk1_to_PhasedXRz);
}

Transform rebase_quil() {
  return rebase_factory(
      {OpType::CZ, OpType::Rx, OpType::Rz}, CircPool::H_CZ_H(),
      CircPool::tk1_to_rzrx);
}

Transform rebase_pyzx() {
  OpTypeSet pyzx_gates = {OpType::SWAP, OpType::CX, OpType::CZ, OpType::H,
                          OpType::X,    OpType::Z,  OpType::S,  OpType::T,
                          OpType::Rx,   OpType::Rz};
  return rebase_factory(pyzx_gates, CircPool::CX(), CircPool::tk1_to_rzrx);
}

Transform rebase_projectq() {
  OpTypeSet projectq_gates = {OpType::SWAP, OpType::CRz, OpType::CX, OpType::CZ,
                              OpType::H,    OpType::X,   OpType::Y,  OpType::Z,
                              OpType::S,    OpType::T,   OpType::V,  OpType::Rx,
                              OpType::Ry,   OpType::Rz};
  return rebase_factory(projectq_gates, CircPool::CX(), CircPool::tk1_to_rzrx);
}

Transform rebase_UFR() {
  return rebase_factory(
      {OpType::CX, OpType::Rz, OpType::H}, CircPool::CX(),
      CircPool::tk1_to_rzh);
}

Transform rebase_OQC() {
  return rebase_factory(
      {OpType::ECR, OpType::Rz, OpType::SX}, CircPool::CX_using_ECR(),
      CircPool::tk1_to_rzsx);
}

// Multiqs: ZZMax
// Singleqs: Rz, PhasedX
Transform rebase_HQS() {
  return rebase_factory(
      {OpType::ZZMax, OpType::Rz, OpType::PhasedX}, CircPool::CX_using_ZZMax(),
      CircPool::tk1_to_PhasedXRz);
}

// Multiqs: TK2
// Singleqs: TK1
Transform rebase_TK() {
  return rebase_factory(
      {OpType::TK2, OpType::TK1}, CircPool::CX_using_TK2(),
      CircPool::tk1_to_tk1);
}

// Multiqs: XXPhase
// Singleqs: Rz, PhasedX
Transform rebase_UMD() {
  return rebase_factory(
      {OpType::XXPhase, OpType::Rz, OpType::PhasedX},
      CircPool::CX_using_XXPhase_0(), CircPool::tk1_to_PhasedXRz);
}

}  // namespace Transforms

}  // namespace tket

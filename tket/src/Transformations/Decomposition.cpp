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

#include "Decomposition.hpp"

#include <functional>
#include <optional>
// replace with c++20 <ranges> when available
#include <boost/range/adaptor/filtered.hpp>

#include "Architecture/Architecture.hpp"
#include "BasicOptimisation.hpp"
#include "Circuit/CircPool.hpp"
#include "Converters/PhasePoly.hpp"
#include "Gate/GatePtr.hpp"
#include "OpType/OpType.hpp"
#include "OpType/OpTypeFunctions.hpp"
#include "Ops/OpPtr.hpp"
#include "OptimisationPass.hpp"
#include "PhasedXFrontier.hpp"
#include "Rebase.hpp"
#include "Replacement.hpp"
#include "Transform.hpp"
#include "Utils/Constants.hpp"
#include "Utils/Expression.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

namespace Transforms {

static bool convert_to_zxz(Circuit &circ);
static bool convert_to_zyz(Circuit &circ);
static bool convert_to_xyx(Circuit &circ);

/**
 * Decompose all multi-qubit unitary gates into TK2 and single-qubit gates.
 *
 * This function does not decompose boxes.
 */
static bool convert_multiqs_TK2(Circuit &circ) {
  bool success = false;
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpType optype = op->get_type();
    if (is_gate_type(optype) && !is_projective_type(optype) &&
        op->n_qubits() >= 2 && (optype != OpType::TK2)) {
      Circuit in_circ = TK2_circ_from_multiq(op);
      Subcircuit sub = {
          {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
      bin.push_back(v);
      circ.substitute(in_circ, sub, Circuit::VertexDeletion::No);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

/**
 * Decompose all multi-qubit unitary gates into CX and single-qubit gates.
 *
 * This function does not decompose boxes.
 */
static bool convert_multiqs_CX(Circuit &circ) {
  bool success = false;
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpType optype = op->get_type();
    if (is_gate_type(optype) && !is_projective_type(optype) &&
        op->n_qubits() >= 2 && (optype != OpType::CX)) {
      Circuit in_circ = CX_circ_from_multiq(op);
      Subcircuit sub = {
          {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
      bin.push_back(v);
      circ.substitute(in_circ, sub, Circuit::VertexDeletion::No);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

static bool convert_to_zxz(Circuit &circ) {
  bool success =
      (decompose_single_qubits_TK1() >> decompose_tk1_to_rzrx()).apply(circ);
  return success;
}

static bool convert_to_zyz(Circuit &circ) {
  static const Expr half = SymEngine::div(Expr(1), Expr(2));
  bool success = decompose_single_qubits_TK1().apply(circ);
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges(v) != 1) continue;
    const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    if (op->get_type() == OpType::TK1) {
      std::vector<Expr> params = op->get_params();
      Circuit replacement(1);
      Expr a = params[2] + half;
      Expr b = params[1];
      Expr c = params[0] - half;
      if (!equiv_0(a, 4)) {
        replacement.add_op<unsigned>(OpType::Rz, a, {0});
      }
      if (!equiv_0(b, 4)) {
        replacement.add_op<unsigned>(OpType::Ry, b, {0});
      }
      if (!equiv_0(c, 4)) {
        replacement.add_op<unsigned>(OpType::Rz, c, {0});
      }
      Subcircuit sub = {
          {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
      bin.push_back(v);
      circ.substitute(replacement, sub, Circuit::VertexDeletion::No);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

static bool convert_to_xyx(Circuit &circ) {
  static const Expr half = SymEngine::div(Expr(1), Expr(2));
  bool success = decompose_single_qubits_TK1().apply(circ);
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges(v) != 1) continue;
    const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    if (op->get_type() == OpType::TK1) {
      std::vector<Expr> params = op->get_params();
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Ry, half, {0});
      replacement.add_op<unsigned>(OpType::Rx, params[2] + half, {0});
      replacement.add_op<unsigned>(OpType::Ry, params[1], {0});
      replacement.add_op<unsigned>(OpType::Rx, params[0] - half, {0});
      replacement.add_op<unsigned>(OpType::Ry, -half, {0});
      remove_redundancies().apply(replacement);
      Subcircuit sub = {
          {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
      bin.push_back(v);
      circ.substitute(replacement, sub, Circuit::VertexDeletion::No);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
  return success;
}

Transform decompose_multi_qubits_TK2() {
  return Transform(convert_multiqs_TK2);
}

Transform decompose_multi_qubits_CX() { return Transform(convert_multiqs_CX); }

static bool convert_singleqs_TK1(Circuit &circ) {
  bool success = false;
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpType optype = op->get_type();
    if (is_gate_type(optype) && !is_projective_type(optype) &&
        op->n_qubits() == 1 && optype != OpType::TK1) {
      std::vector<Expr> tk1_angs = as_gate_ptr(op)->get_tk1_angles();
      Circuit rep(1);
      rep.add_op<unsigned>(
          OpType::TK1, {tk1_angs[0], tk1_angs[1], tk1_angs[2]}, {0});
      circ.substitute(rep, v, Circuit::VertexDeletion::No);
      circ.add_phase(tk1_angs[3]);
      bin.push_back(v);
      success = true;
    }
  }
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::Yes);
  return success;
}

Transform decompose_single_qubits_TK1() {
  return Transform(convert_singleqs_TK1);
}

Transform decompose_ZYZ_to_TK1() {
  return Transform([](Circuit &circ) {
    bool success = false;
    static const Expr zero(0);
    static const Expr half = SymEngine::div(Expr(1), Expr(2));
    VertexList bin;
    VertexVec inputs = circ.q_inputs();
    for (VertexVec::iterator i = inputs.begin(); i != inputs.end(); ++i) {
      Edge e = circ.get_nth_out_edge(*i, 0);
      Vertex v = circ.target(e);
      while (!is_final_q_type(circ.get_OpType_from_Vertex(v))) {
        if (circ.get_OpType_from_Vertex(v) == OpType::Rz) {
          const Op_ptr v_g = circ.get_Op_ptr_from_Vertex(v);
          Expr angle_1 = v_g->get_params()[0];
          Edge e1 = circ.get_next_edge(v, e);
          Vertex v2 = circ.target(e1);
          if (circ.get_OpType_from_Vertex(v2) == OpType::Ry) {
            const Op_ptr v2_g = circ.get_Op_ptr_from_Vertex(v2);
            Expr angle_2 = v2_g->get_params()[0];
            Edge e2 = circ.get_next_edge(v2, e1);
            Vertex v3 = circ.target(e2);
            bin.push_back(v2);
            circ.remove_vertex(
                v2, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            Expr angle_3 = zero;
            if (circ.get_OpType_from_Vertex(v3) == OpType::Rz) {
              const Op_ptr v3_g = circ.get_Op_ptr_from_Vertex(v3);
              angle_3 = v3_g->get_params()[0];
              circ.remove_vertex(
                  v3, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
              bin.push_back(v3);
            }
            std::vector<Expr> new_params = {
                angle_3 + half, angle_2, angle_1 - half};
            circ.dag[v] = {get_op_ptr(OpType::TK1, new_params)};
          } else {
            circ.dag[v] = {get_op_ptr(OpType::TK1, {zero, zero, angle_1})};
          }
        } else if (circ.get_OpType_from_Vertex(v) == OpType::Ry) {
          const Op_ptr v_g = circ.get_Op_ptr_from_Vertex(v);
          Expr angle_2 = v_g->get_params()[0];
          Expr angle_3 = zero;
          Edge e1 = circ.get_next_edge(v, e);
          Vertex v2 = circ.target(e1);
          if (circ.get_OpType_from_Vertex(v2) == OpType::Rz) {
            const Op_ptr v2_g = circ.get_Op_ptr_from_Vertex(v2);
            angle_3 = v2_g->get_params()[0];
            circ.remove_vertex(
                v2, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            bin.push_back(v2);
          }
          std::vector<Expr> new_params = {angle_3 + half, angle_2, -half};
          circ.dag[v] = {get_op_ptr(OpType::TK1, new_params)};
        }
        e = circ.get_next_edge(v, e);
        v = circ.target(e);
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform decompose_ZX() { return Transform(convert_to_zxz); }

Transform decompose_ZY() { return Transform(convert_to_zyz); }

Transform decompose_XY() { return Transform(convert_to_xyx); }

Transform decompose_tk1_to_rzrx() {
  return Transform([](Circuit &circ) {
    bool success = false;
    auto [it, end] = boost::vertices(circ.dag);
    for (auto next = it; it != end; it = next) {
      ++next;
      if (circ.get_OpType_from_Vertex(*it) == OpType::TK1) {
        success = true;
        const Op_ptr g = circ.get_Op_ptr_from_Vertex(*it);
        const std::vector<Expr> &params = g->get_params();
        Circuit newcirc =
            CircPool::tk1_to_rzrx(params[0], params[1], params[2]);
        Subcircuit sc = {
            {circ.get_in_edges(*it)}, {circ.get_all_out_edges(*it)}, {*it}};
        circ.substitute(newcirc, sc, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform decompose_CX_to_ECR() {
  return Transform([](Circuit &circ) {
    bool success = false;
    auto [i, end] = boost::vertices(circ.dag);
    for (auto next = i; i != end; i = next) {
      ++next;
      if (circ.get_OpType_from_Vertex(*i) == OpType::CX) {
        success = true;
        Subcircuit sub{circ.get_in_edges(*i), circ.get_all_out_edges(*i), {*i}};
        circ.substitute(
            CircPool::CX_using_ECR(), sub, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform decompose_CX_to_HQS2() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::CX) {
        success = true;
        bin.push_back(v);
        Subcircuit sub = {circ.get_in_edges(v), circ.get_all_out_edges(v)};
        circ.substitute(
            CircPool::CX_using_ZZMax(), sub, Circuit::VertexDeletion::No);
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

/* --Rz(a)--Rx(b)--R(c)-- => --Rz(a+c)--PhasedX(b,c)-- */
Transform decompose_ZX_to_HQS1() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList to_bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::Rx) {
        success = true;
        const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
        Expr theta = g->get_params()[0];
        Vertex prev_vert = *circ.get_predecessors(v).begin();
        Vertex next_vert = *circ.get_successors(v).begin();
        if (circ.get_OpType_from_Vertex(prev_vert) == OpType::Rz &&
            circ.get_OpType_from_Vertex(next_vert) == OpType::Rz) {
          const Op_ptr prev_g = circ.get_Op_ptr_from_Vertex(prev_vert);
          const Op_ptr next_g = circ.get_Op_ptr_from_Vertex(next_vert);
          Expr phi = next_g->get_params()[0];
          std::vector<Expr> params{theta, phi};
          circ.dag[v].op = get_op_ptr(OpType::PhasedX, params);
          circ.remove_vertex(
              next_vert, Circuit::GraphRewiring::Yes,
              Circuit::VertexDeletion::No);
          to_bin.push_back(next_vert);
          Expr new_param = prev_g->get_params()[0] + phi;
          circ.dag[prev_vert].op = get_op_ptr(OpType::Rz, new_param);
        } else {
          // if no Rz, initialise a PhasedX op with theta=Rx.params[0],phi=0
          Expr phi(0);
          std::vector<Expr> params{theta, phi};
          circ.dag[v].op = get_op_ptr(OpType::PhasedX, params);
        }
      }
    }
    circ.remove_vertices(
        to_bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    remove_redundancies().apply(circ);
    return success;
  });
}

// Decompose CX into MolmerSorensen as:
// ---C---         -V-S-|-H-
//    |      -->    XX(pi/4)
// ---X---         -----|-Vdg-
Transform decompose_MolmerSorensen() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) != OpType::CX) continue;
      EdgeVec outs = circ.get_all_out_edges(v);
      if (outs.size() == 2) {
        Vertex next = circ.target(outs[0]);
        // Is the next operation equivalent to an Rx, up to phase?
        Op_ptr next_g = circ.get_Op_ptr_from_Vertex(next);
        OpType next_type = next_g->get_type();
        if (is_single_qubit_type(next_type) && !is_projective_type(next_type)) {
          std::vector<Expr> angles = as_gate_ptr(next_g)->get_tk1_angles();
          if (equiv_0(angles[0], 2) && equiv_0(angles[2], 2)) {
            Expr angle = angles[1];
            Expr phase = angles[3];
            if (!equiv_0(angles[0], 4)) phase += 1;
            if (!equiv_0(angles[2], 4)) phase += 1;
            Edge next_e = circ.get_nth_out_edge(next, 0);
            Vertex last = circ.target(next_e);
            if (circ.get_OpType_from_Vertex(last) == OpType::CX &&
                circ.get_nth_in_edge(last, 1) == outs[1]) {
              // Recognise exp(-i XX * angle * pi/2)
              const Op_ptr op_ptr = get_op_ptr(OpType::XXPhase, angle);
              circ.dag[v] = {op_ptr};
              bin.push_back(next);
              circ.remove_vertex(
                  next, Circuit::GraphRewiring::Yes,
                  Circuit::VertexDeletion::No);
              bin.push_back(last);
              circ.remove_vertex(
                  last, Circuit::GraphRewiring::Yes,
                  Circuit::VertexDeletion::No);
              circ.add_phase(phase);
              success = true;
              continue;
            }
          }
        }
        // Replace remaining CX gates
        Subcircuit sub = {{circ.get_in_edges(v)}, {outs}, {v}};
        bin.push_back(v);
        circ.substitute(
            CircPool::CX_using_XXPhase_1(), sub, Circuit::VertexDeletion::No);
        success = true;
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

static double get_ZZPhase_fidelity(
    const std::array<double, 3> &k, unsigned nb_cx) {
  switch (nb_cx) {
    case 0:
      return trace_fidelity(k[0], k[1], k[2]);
    case 1:
      return trace_fidelity(0, k[1], k[2]);
    case 2:
      return trace_fidelity(0, 0, k[2]);
    default:
      return 1.;
  }
}

// Try to decompose a TK2 gate using different gate sets, find the one with
// the highest fidelity.
// If no fidelities are provided, (best_optype, n_gates) is left unchanged.
static void best_noise_aware_decomposition(
    const std::array<double, 3> &angles, const TwoQbFidelities &fid,
    OpType &best_optype, unsigned &n_gates) {
  double max_fid = 0.;

  // Try decomposition using CX or equivalent gates.
  double cx_fid = std::max(
      fid.CX_fidelity ? fid.CX_fidelity.value() : 0.,
      fid.ZZMax_fidelity ? fid.ZZMax_fidelity.value() : 0.);
  bool zzmax_is_better = false;
  if (cx_fid < EPS) {
    if (!fid.ZZPhase_fidelity) {
      // No fidelity is defined, so default to CX
      cx_fid = 1.;
    }
  } else {
    zzmax_is_better = fid.CX_fidelity < fid.ZZMax_fidelity;
  }
  if (cx_fid > EPS) {
    for (unsigned n_cx = 0; n_cx <= 3; ++n_cx) {
      double ncx_fid = get_CX_fidelity(angles, n_cx) * pow(cx_fid, n_cx);
      if (ncx_fid > max_fid) {
        max_fid = ncx_fid;
        best_optype = zzmax_is_better ? OpType::ZZMax : OpType::CX;
        n_gates = n_cx;
      }
    }
  }

  // Try decomposition using ZZPhase(α)
  if (fid.ZZPhase_fidelity) {
    double zz_fid = 1.;
    // If ZZMax is available, ZZPhase is only interesting when used once.
    // (two ZZPhase can always be written using two ZZmax)
    unsigned max_nzz = fid.ZZMax_fidelity ? 1 : 3;
    for (unsigned n_zz = 0; n_zz <= max_nzz; ++n_zz) {
      if (n_zz > 0) {
        double gate_fid = (*fid.ZZPhase_fidelity)(angles[n_zz - 1]);
        if (gate_fid < 0 || gate_fid > 1) {
          throw NotValid(
              "ZZPhase_fidelity returned a value outside of [0, 1].");
        }
        zz_fid *= gate_fid;
      }
      double nzz_fid = get_ZZPhase_fidelity(angles, n_zz) * zz_fid;
      if (nzz_fid > max_fid) {
        max_fid = nzz_fid;
        best_optype = OpType::ZZPhase;
        n_gates = n_zz;
      }
    }
  }
}

// Try to decompose a TK2 gate exactly using different gate sets.
// The fidelities are used as an indication of which gate set is usable, but
// the actual values of the fidelities will be ignored.
//
// Relies on default values of best_optype and n_gates if no optimisation can be
// performed.
static void best_exact_decomposition(
    const std::array<Expr, 3> &angles, const TwoQbFidelities &fid,
    OpType &best_optype, unsigned &n_gates) {
  // Prefer CX/ZZMax when possible.
  if (fid.CX_fidelity || fid.ZZMax_fidelity) {
    bool zzmax_is_better =
        !fid.CX_fidelity || fid.CX_fidelity < fid.ZZMax_fidelity;
    best_optype = zzmax_is_better ? OpType::ZZMax : OpType::CX;
  } else if (fid.ZZPhase_fidelity) {
    // Only ZZPhase has fidelity, so use ZZPhase.
    best_optype = OpType::ZZPhase;
  }

  if (best_optype == OpType::CX || best_optype == OpType::ZZMax) {
    // Reduce n_gates if possible.
    if (equiv_0(angles[2], 4)) {
      n_gates = 2;
    }
  } else if (best_optype == OpType::ZZPhase) {
    // Reduce n_gates if possible.
    if (equiv_0(angles[2], 4)) {
      n_gates = 2;
      if (equiv_0(angles[1], 4)) {
        n_gates = 1;
      }
    }
  }

  // Finally, handle the only case where ZZPhase is preferable over ZZMax.
  if (fid.ZZPhase_fidelity && equiv_0(angles[2], 4) && equiv_0(angles[1], 4) &&
      n_gates > 1) {
    n_gates = 1;
    best_optype = OpType::ZZPhase;
  }
}

// Whether the TK2 angles are normalised.
//
// Numerical values must be in the Weyl chamber, ie 1/2 >= k_x >= k_y >= |k_z|.
// Symbolic values must come before any numerical value in the array.
static bool in_weyl_chamber(const std::array<Expr, 3> &k) {
  bool is_symbolic = true;
  double last_val = .5;
  for (unsigned i = 0; i < k.size(); ++i) {
    std::optional<double> eval = eval_expr_mod(k[i], 4);
    if (eval) {
      is_symbolic = false;
      if (i + 1 == k.size()) {
        double abs_eval = std::min(*eval, -(*eval) + 4);
        if (abs_eval > last_val) {
          return false;
        }
      } else {
        if (*eval > last_val) {
          return false;
        }
      }
      last_val = *eval;
    } else if (!is_symbolic) {
      return false;
    }
  }
  return true;
}

/**
 * @brief TK2 expressed (approximately) as CX/ZZMax or ZZPhase.
 *
 * This is the core logic of how to decompose a TK2 gate into other two-qubit
 * gates, taking hardware fidelities into account for optimal approximate
 * decompositions.
 *
 * Decomposes to whatever gate type has non-nullopt fidelity. If there are
 * multiple options, choose the best one. Defaults to CX if no fidelities are
 * provided.
 *
 * Symbolic parameters are supported. In that case, decompositions are exact.
 *
 * @param angles The TK2 parameters
 * @param fid The two-qubit gate fidelities
 * @return Circuit TK2-equivalent circuit
 */
static Circuit TK2_replacement(
    const std::array<Expr, 3> &angles, const TwoQbFidelities &fid) {
  if (!in_weyl_chamber(angles)) {
    throw NotValid("TK2 params are not normalised to Weyl chamber.");
  }
  OpType best_optype = OpType::CX;  // default to using CX
  unsigned n_gates = 3;             // default to 3x CX

  // Try to evaluate exprs to doubles.
  std::array<double, 3> angles_eval;
  unsigned last_angle = 0;
  for (const Expr &e : angles) {
    std::optional<double> eval = eval_expr_mod(e);
    if (eval) {
      angles_eval[last_angle++] = *eval;
    } else {
      break;
    }
  }

  if (last_angle <= 2) {
    // Not all angles could be resolved numerically.
    // For symbolic angles, we can only provide an exact decomposition.
    best_exact_decomposition(angles, fid, best_optype, n_gates);
  } else {
    // For non-symbolic angles, we can find the optimal number of gates
    // using the gate fidelities provided.
    best_noise_aware_decomposition(angles_eval, fid, best_optype, n_gates);
  }

  // Build circuit for substitution.
  Circuit sub(2);
  switch (best_optype) {
    case OpType::ZZMax:
    case OpType::CX: {
      switch (n_gates) {
        case 0:
          break;
        case 1: {
          sub.append(CircPool::approx_TK2_using_1xCX());
          break;
        }
        case 2: {
          sub.append(CircPool::approx_TK2_using_2xCX(angles[0], angles[1]));
          break;
        }
        case 3: {
          sub.append(CircPool::TK2_using_3xCX(angles[0], angles[1], angles[2]));
          break;
        }
        default:
          throw NotValid("Number of CX invalid in decompose_TK2");
      }
      if (best_optype == OpType::ZZMax) {
        decompose_CX_to_HQS2().apply(sub);
      }
      break;
    }
    case OpType::ZZPhase: {
      switch (n_gates) {
        case 0:
          break;
        case 1: {
          sub.append(CircPool::approx_TK2_using_1xZZPhase(angles[0]));
          break;
        }
        case 2: {
          sub.append(
              CircPool::approx_TK2_using_2xZZPhase(angles[0], angles[1]));
          break;
        }
        case 3: {
          sub.append(
              CircPool::TK2_using_ZZPhase(angles[0], angles[1], angles[2]));
          break;
        }
        default:
          throw NotValid("Number of ZZPhase invalid in decompose_TK2");
      }
      break;
    }
    default:
      throw NotValid("Unrecognised target OpType in decompose_TK2");
  }
  return sub;
}

Transform decompose_TK2() { return decompose_TK2({}); }
Transform decompose_TK2(const TwoQbFidelities &fid) {
  if (fid.ZZMax_fidelity) {
    if (*fid.ZZMax_fidelity < 0 || *fid.ZZMax_fidelity > 1) {
      throw NotValid("ZZMax fidelity must be between 0 and 1.");
    }
  }
  if (fid.CX_fidelity) {
    if (*fid.CX_fidelity < 0 || *fid.CX_fidelity > 1) {
      throw NotValid("CX fidelity must be between 0 and 1.");
    }
  }
  if (fid.ZZMax_fidelity && fid.ZZPhase_fidelity) {
    if (*fid.ZZMax_fidelity < (*fid.ZZPhase_fidelity)(.5)) {
      throw NotValid(
          "The ZZMax fidelity cannot be smaller than the ZZPhase(0.5) "
          "fidelity");
    }
  }
  return Transform([fid](Circuit &circ) {
    bool success = false;

    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) != OpType::TK2) continue;

      success = true;
      auto params = circ.get_Op_ptr_from_Vertex(v)->get_params();
      TKET_ASSERT(params.size() == 3);
      std::array<Expr, 3> angles{params[0], params[1], params[2]};

      Circuit sub = TK2_replacement(angles, fid);
      bin.push_back(v);
      circ.substitute(sub, v, Circuit::VertexDeletion::No);
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

    return success;
  });
}

Transform decompose_ZZPhase() {
  return Transform([](Circuit &circ) {
    bool success = decompose_PhaseGadgets().apply(circ);
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      switch (circ.get_OpType_from_Vertex(v)) {
        case OpType::PhaseGadget: {
          success = true;
          const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
          TKET_ASSERT(g->get_params().size() == 1);
          circ.dag[v] = {get_op_ptr(OpType::ZZPhase, g->get_params()[0])};
          break;
        }
        case OpType::XXPhase: {
          success = true;
          const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
          TKET_ASSERT(g->get_params().size() == 1);
          Circuit sub = CircPool::XXPhase_using_ZZPhase(g->get_params()[0]);
          circ.substitute(sub, v, Circuit::VertexDeletion::No);
          bin.push_back(v);
          break;
        }
        case OpType::YYPhase: {
          success = true;
          const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
          TKET_ASSERT(g->get_params().size() == 1);
          Circuit sub = CircPool::YYPhase_using_ZZPhase(g->get_params()[0]);
          circ.substitute(sub, v, Circuit::VertexDeletion::No);
          bin.push_back(v);
          break;
        }
        default:
          break;
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

/**
 * Specification of a sequence of Clifford gates and a phase
 *
 * The order of the gates is (Z)(X)(S)(V)(S). Each int is 1 or 0 representing
 * the presence or absence of the corresponding gate.
 */
typedef struct {
  int Z0;
  int X1;
  int S2;
  int V3;
  int S4;
  double p;
} std_cliff_spec_t;

/**
 * The (i,j,k) entry in this table represents TK1(i/2, j/2, k/2).
 *
 * Where there is more than one decomposition the number of gates is minimized.
 */
static const std_cliff_spec_t tk1_table[4][4][4] = {
    {
        {
            {0, 0, 0, 0, 0, 0.0},
            {0, 0, 0, 0, 1, -0.25},
            {1, 0, 0, 0, 0, -0.5},
            {1, 0, 0, 0, 1, -0.75},
        },
        {
            {0, 0, 0, 1, 0, 0.0},
            {0, 0, 1, 1, 0, -0.25},
            {1, 0, 0, 1, 0, -0.5},
            {1, 0, 1, 1, 0, -0.75},
        },
        {
            {0, 1, 0, 0, 0, -0.5},
            {1, 1, 0, 0, 1, 0.75},
            {1, 1, 0, 0, 0, 1.0},
            {0, 1, 0, 0, 1, 0.25},
        },
        {
            {0, 1, 0, 1, 0, -0.5},
            {1, 1, 1, 1, 0, 0.75},
            {1, 1, 0, 1, 0, 1.0},
            {0, 1, 1, 1, 0, 0.25},
        },
    },
    {
        {
            {0, 0, 0, 0, 1, -0.25},
            {1, 0, 0, 0, 0, -0.5},
            {1, 0, 0, 0, 1, -0.75},
            {0, 0, 0, 0, 0, -1.0},
        },
        {
            {0, 0, 0, 1, 1, -0.25},
            {0, 0, 1, 1, 1, -0.5},
            {1, 0, 0, 1, 1, -0.75},
            {1, 0, 1, 1, 1, -1.0},
        },
        {
            {0, 1, 0, 0, 1, -0.75},
            {0, 1, 0, 0, 0, -0.5},
            {1, 1, 0, 0, 1, 0.75},
            {1, 1, 0, 0, 0, 1.0},
        },
        {
            {0, 1, 0, 1, 1, -0.75},
            {1, 1, 1, 1, 1, 0.5},
            {1, 1, 0, 1, 1, 0.75},
            {0, 1, 1, 1, 1, 0.0},
        },
    },
    {
        {
            {1, 0, 0, 0, 0, -0.5},
            {1, 0, 0, 0, 1, -0.75},
            {0, 0, 0, 0, 0, -1.0},
            {0, 0, 0, 0, 1, 0.75},
        },
        {
            {1, 1, 0, 1, 0, 0.0},
            {0, 1, 1, 1, 0, -0.75},
            {0, 1, 0, 1, 0, -0.5},
            {1, 1, 1, 1, 0, 0.75},
        },
        {
            {1, 1, 0, 0, 0, 0.0},
            {0, 1, 0, 0, 1, -0.75},
            {0, 1, 0, 0, 0, -0.5},
            {1, 1, 0, 0, 1, 0.75},
        },
        {
            {1, 0, 0, 1, 0, 0.5},
            {1, 0, 1, 1, 0, 0.25},
            {0, 0, 0, 1, 0, 0.0},
            {0, 0, 1, 1, 0, -0.25},
        },
    },
    {
        {
            {1, 0, 0, 0, 1, -0.75},
            {0, 0, 0, 0, 0, -1.0},
            {0, 0, 0, 0, 1, 0.75},
            {1, 0, 0, 0, 0, 0.5},
        },
        {
            {1, 1, 0, 1, 1, -0.25},
            {0, 1, 1, 1, 1, -1.0},
            {0, 1, 0, 1, 1, -0.75},
            {1, 1, 1, 1, 1, 0.5},
        },
        {
            {1, 1, 0, 0, 1, -0.25},
            {1, 1, 0, 0, 0, 0.0},
            {0, 1, 0, 0, 1, -0.75},
            {0, 1, 0, 0, 0, -0.5},
        },
        {
            {1, 0, 0, 1, 1, 0.25},
            {1, 0, 1, 1, 1, 0.0},
            {0, 0, 0, 1, 1, -0.25},
            {0, 0, 1, 1, 1, -0.5},
        },
    },
};

/**
 * Clifford circuit equivalent to TK1(i/2, j/2, k/2)
 *
 * @pre 0 <= i, j, k < 8
 * @post circuit consists of V, S, X and Z gates only, in order (Z)(X)(S)(V)(S)
 */
static Circuit clifford_from_tk1(int i, int j, int k) {
  std_cliff_spec_t spec = tk1_table[i % 4][j % 4][k % 4];
  if (i >= 4) spec.p += 1;
  if (j >= 4) spec.p += 1;
  if (k >= 4) spec.p += 1;

  Circuit c(1);
  if (spec.Z0) {
    c.add_op<unsigned>(OpType::Z, {0});
  }
  if (spec.X1) {
    c.add_op<unsigned>(OpType::X, {0});
  }
  if (spec.S2) {
    c.add_op<unsigned>(OpType::S, {0});
  }
  if (spec.V3) {
    c.add_op<unsigned>(OpType::V, {0});
  }
  if (spec.S4) {
    c.add_op<unsigned>(OpType::S, {0});
  }
  c.add_phase(spec.p);

  return c;
}

Transform decompose_cliffords_std() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType type = op->get_type();
      if (type != OpType::V && type != OpType::S && type != OpType::X &&
          type != OpType::Z && is_single_qubit_unitary_type(type) &&
          op->is_clifford()) {
        std::vector<Expr> tk1_param_exprs = as_gate_ptr(op)->get_tk1_angles();
        bool all_reduced = true;
        bool all_roundable = true;
        std::vector<int> iangles(3);
        for (int i = 0; i < 3; i++) {
          std::optional<double> reduced = eval_expr_mod(tk1_param_exprs[i], 4);
          if (!reduced)
            all_reduced = false;
          else {
            double angle = 2 * reduced.value();  // > 0
            iangles[i] = int(angle + 0.5);       // nearest integer
            if (std::abs(angle - iangles[i]) >= EPS) {
              all_roundable = false;
            }
            iangles[i] %= 8;  // 8 --> 0
          }
        }
        if (!(all_reduced && all_roundable)) continue;
        Circuit replacement =
            clifford_from_tk1(iangles[0], iangles[1], iangles[2]);
        Subcircuit sub = {
            {circ.get_in_edges(v)}, {circ.get_all_out_edges(v)}, {v}};
        bin.push_back(v);
        circ.substitute(replacement, sub, Circuit::VertexDeletion::No);
        circ.add_phase(tk1_param_exprs[3]);
        success = true;
      } else if (type == OpType::TK2 && op->is_clifford()) {
        auto params = op->get_params();
        TKET_ASSERT(params.size() == 3);
        // TODO: Maybe handle TK2 gates natively within clifford_simp?
        Circuit replacement =
            CircPool::TK2_using_CX(params[0], params[1], params[2]);
        decompose_cliffords_std().apply(replacement);
        bin.push_back(v);
        circ.substitute(replacement, v, Circuit::VertexDeletion::No);
        success = true;
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform decompose_ZX_to_cliffords() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      const Op_ptr op_ptr = circ.get_Op_ptr_from_Vertex(v);
      if (op_ptr->get_type() == OpType::Rz ||
          op_ptr->get_type() == OpType::Rx) {
        Expr param = op_ptr->get_params()[0];
        std::optional<double> reduced = eval_expr_mod(param, 4);
        bool roundable = false;
        int iangle = 0;
        if (reduced) {
          double angle = 2 * reduced.value();  // >= 0
          iangle = int(angle + 0.5);           // nearest integer
          iangle %= 8;                         // {0,1,2,3,4,5,6,7}
          roundable = (std::abs(angle - iangle) < EPS);
        }
        if (roundable) {
          bool is_rz = op_ptr->get_type() == OpType::Rz;
          switch (iangle % 4) {
            case 0: {
              bin.push_back(v);
              circ.remove_vertex(
                  v, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
              break;
            }
            case 1: {
              if (is_rz) {
                circ.dag[v] = {get_op_ptr(OpType::S)};
                circ.add_phase(-0.25);
              } else {
                circ.dag[v] = {get_op_ptr(OpType::V)};
              }
              break;
            }
            case 2: {
              if (is_rz) {
                circ.dag[v] = {get_op_ptr(OpType::Z)};
              } else {
                circ.dag[v] = {get_op_ptr(OpType::X)};
              }
              circ.add_phase(-0.5);
              break;
            }
            case 3: {
              if (is_rz) {
                circ.dag[v] = {get_op_ptr(OpType::Sdg)};
                circ.add_phase(-0.75);
              } else {
                circ.dag[v] = {get_op_ptr(OpType::Vdg)};
                circ.add_phase(1);
              }
              break;
            }
          }
          if (iangle >= 4) circ.add_phase(1);
          success = true;
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform decompose_PhaseGadgets() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList big_bin;
    auto [it, end] = boost::vertices(circ.dag);
    for (auto next = it; it != end; it = next) {
      ++next;
      if (circ.get_OpType_from_Vertex(*it) == OpType::CX &&
          circ.n_out_edges(*it) == 2) {
        EdgeVec outs = circ.get_all_out_edges(*it);
        Vertex next_v = circ.target(outs[1]);
        Op_ptr g = circ.get_Op_ptr_from_Vertex(next_v);
        OpType type = g->get_type();
        if (type == OpType::Rz || type == OpType::U1 ||
            (type == OpType::TK1 && equiv_0(g->get_params()[1]))) {
          Vertex last_v = circ.get_next_pair(next_v, outs[1]).first;
          if (circ.get_OpType_from_Vertex(last_v) == OpType::CX &&
              circ.get_nth_in_edge(last_v, 0) == outs[0]) {
            VertexList bin = {next_v, last_v};
            big_bin.push_back(next_v);
            big_bin.push_back(last_v);
            circ.remove_vertices(
                bin, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            Expr t = g->get_params()[0];
            if (type == OpType::TK1) {
              t += g->get_params()[2];
            }
            circ.dag[*it] = {get_op_ptr(OpType::PhaseGadget, {t}, 2)};
            if (type == OpType::U1) {
              circ.add_phase(t / 2);
            } else if (
                type == OpType::TK1 && equiv_val(g->get_params()[1], 2, 4)) {
              circ.add_phase(1);
            }
            success = true;
          }
        }
        if (type == OpType::CX) {
          if (circ.get_target_port(outs[1]) == 1) {
            Vertex rx = circ.source(circ.get_nth_in_edge(next_v, 0));
            if (circ.get_OpType_from_Vertex(rx) == OpType::Rx) {
              if (rx == circ.target(outs[0])) {
                const Op_ptr rx_g = (circ.get_Op_ptr_from_Vertex(rx));
                VertexList bin = {rx, next_v};
                big_bin.push_back(next_v);
                big_bin.push_back(rx);
                Circuit replacement(2);
                circ.remove_vertices(
                    bin, Circuit::GraphRewiring::Yes,
                    Circuit::VertexDeletion::No);
                replacement.add_op<unsigned>(OpType::H, {0});
                replacement.add_op<unsigned>(OpType::H, {1});
                replacement.add_op<unsigned>(
                    OpType::PhaseGadget, rx_g->get_params()[0], {0, 1});
                replacement.add_op<unsigned>(OpType::H, {0});
                replacement.add_op<unsigned>(OpType::H, {1});
                Subcircuit sub = {
                    circ.get_in_edges(*it), circ.get_all_out_edges(*it), {*it}};
                circ.substitute(replacement, sub, Circuit::VertexDeletion::Yes);
                success = true;
              }
            }
          }
        }
      }
    }
    circ.remove_vertices(
        big_bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform decomp_boxes() {
  return Transform([](Circuit &circ) { return circ.decompose_boxes(); });
}

Transform compose_phase_poly_boxes(const unsigned min_size) {
  return Transform([=](Circuit &circ) {
    // replace wireswaps with three CX
    while (circ.has_implicit_wireswaps()) {
      qubit_map_t perm = circ.implicit_qubit_permutation();
      for (const std::pair<const Qubit, Qubit> &pair : perm) {
        if (pair.first != pair.second) {
          circ.replace_implicit_wire_swap(pair.first, pair.second);
          break;
        }
      }
    }

    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ, min_size);
    conv.convert();
    circ = conv.get_circuit();
    return true;
  });
}

Transform decompose_SWAP(const Circuit &replacement_circuit) {
  return Transform([=](Circuit &circ) {
    if (!replacement_circuit.is_simple()) throw SimpleOnly();
    return circ.substitute_all(replacement_circuit, get_op_ptr(OpType::SWAP));
  });
}

static void swap_sub(
    Circuit &circ, const Circuit &swap_circ_1, const Circuit &swap_circ_2,
    Subcircuit &sub, const std::pair<port_t, port_t> &port_comp) {
  std::pair<port_t, port_t> comp = {0, 1};
  // Ports only come in 2 cases, {0,1} or {1,0}. if {0,1} (first case),
  // swap_circ_1 leaves a CX{0,1} next to current CX{0,1}, if not we can
  // assume second case.
  if (port_comp == comp)
    circ.substitute(swap_circ_1, sub, Circuit::VertexDeletion::Yes);
  else
    circ.substitute(swap_circ_2, sub, Circuit::VertexDeletion::Yes);
}

Transform decompose_SWAP_to_CX(const Architecture &arc) {
  // Note that the default argument will be out of scope at call-time!
  //  => we replace the default empty Architecture with nullptr
  // we need to keep arc as a pointer as there is no such thing as
  // optional references in std
  const Architecture *arc_ptr = arc.n_nodes() ? &arc : nullptr;
  return Transform([arc_ptr](Circuit &circ) {
    bool success = false;
    std::vector<std::pair<Vertex, bool>> bin;
    for (Circuit::CommandIterator it = circ.begin(); it != circ.end(); ++it) {
      if (it->get_op_ptr()->get_type() == OpType::SWAP) {
        unit_vector_t qbs = it->get_args();
        node_vector_t nodes = {qbs.begin(), qbs.end()};
        if (arc_ptr != nullptr && arc_ptr->node_exists(nodes[0]) &&
            arc_ptr->node_exists(nodes[1]) &&
            arc_ptr->edge_exists(nodes[1], nodes[0])) {
          bin.push_back({it.get_vertex(), true});
        } else {
          bin.push_back({it.get_vertex(), false});
        }
      }
    }

    for (std::pair<Vertex, bool> v : bin) {
      success = true;
      // Get predecessor vertices and successor vertices and find subcircuit
      // for replacement
      VertexVec preds = circ.get_predecessors(v.first);
      VertexVec succs = circ.get_successors(v.first);
      EdgeVec in_edges = circ.get_in_edges(v.first);
      EdgeVec out_edges = circ.get_all_out_edges(v.first);
      Subcircuit sub = {in_edges, out_edges, {v.first}};

      if (preds.size() == 1 &&
          circ.get_OpType_from_Vertex(preds[0]) == OpType::CX) {
        // if conditions requires that there is a CX gate before SWAP
        swap_sub(
            circ, CircPool::SWAP_using_CX_0(), CircPool::SWAP_using_CX_1(), sub,
            {circ.get_source_port(in_edges[0]),
             circ.get_source_port(in_edges[1])});
      } else if (
          succs.size() == 1 &&
          circ.get_OpType_from_Vertex(succs[0]) == OpType::CX) {
        // if no CX gate before SWAP, check after.
        swap_sub(
            circ, CircPool::SWAP_using_CX_0(), CircPool::SWAP_using_CX_1(), sub,
            {circ.get_target_port(out_edges[0]),
             circ.get_target_port(out_edges[1])});
      } else if (v.second) {
        // We assume in general that a CX gate saving is desired over H gate
        // savings. If SWAP doesn't lend itself to annihlation though, the
        // SWAP is inserted to reduce number of H gates added in a 'directed'
        // CX decomposition. SWAP_using_CX_1 is added if the backwards
        // direction is available on the architecture
        circ.substitute(
            CircPool::SWAP_using_CX_1(), sub, Circuit::VertexDeletion::Yes);
      } else {
        // SWAP_using_CX_0 is added as a default option
        circ.substitute(
            CircPool::SWAP_using_CX_0(), sub, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform decompose_BRIDGE_to_CX() {
  return Transform([](Circuit &circ) {
    bool success = false;
    // Collect all BRIDGE type vertices
    std::vector<std::pair<Vertex, bool>> bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::BRIDGE) {
        bin.push_back({v, false});
      }
      if (circ.get_OpType_from_Vertex(v) == OpType::Conditional) {
        const Conditional &b =
            static_cast<const Conditional &>(*circ.get_Op_ptr_from_Vertex(v));
        if (b.get_op()->get_type() == OpType::BRIDGE) {
          bin.push_back({v, true});
        }
      }
    }

    auto BRIDGE_sub =
        [&circ](std::pair<Vertex, bool> &candidate, Circuit BRIDGE_circ) {
          if (candidate.second)
            circ.substitute_conditional(
                BRIDGE_circ, candidate.first, Circuit::VertexDeletion::Yes);
          else
            circ.substitute(
                BRIDGE_circ, candidate.first, Circuit::VertexDeletion::Yes);
        };

    for (std::pair<Vertex, bool> v : bin) {
      // Get predecessor vertices and successor vertices and find subcircuit
      // for replacement
      success = true;
      VertexVec preds = circ.get_predecessors(v.first);
      VertexVec succs = circ.get_successors(v.first);
      EdgeVec in_edges = circ.get_in_edges(v.first);
      EdgeVec out_edges = circ.get_all_out_edges(v.first);
      Subcircuit sub = {in_edges, out_edges, {v.first}};

      bool done = false;
      if (preds.size() < 3) {  // Implies some predecessor vertices are in a
                               // multi-qubit op together
        VertexVec comps = {
            circ.source(in_edges[0]), circ.source(in_edges[1]),
            circ.source(in_edges[2])};
        if (comps[0] ==
            comps[1]) {  // First two qubits in BRIDGE in multi-qubit op
          BRIDGE_sub(v, CircPool::BRIDGE_using_CX_0());
          done = true;
        } else if (comps[2] == comps[1]) {  // Second two qubits in BRIDGE in
                                            // multi-qubit op together before
                                            // BRIDGE
          BRIDGE_sub(v, CircPool::BRIDGE_using_CX_1());
          done = true;
        }
      }
      if (succs.size() < 3 && !done) {  // Implies some successor vertices are
                                        // in a multiq-ubit op together
        VertexVec comps = {
            circ.target(out_edges[0]), circ.target(out_edges[1]),
            circ.target(out_edges[2])};
        if (comps[0] == comps[1]) {  // First two qubits in BRIDGE in
                                     // multi-qubit op together before BRIDGE
          BRIDGE_sub(v, CircPool::BRIDGE_using_CX_1());
          done = true;
        } else if (comps[2] == comps[1]) {  // Second two qubits in BRIDGE in
                                            // multi-qubit op together before
                                            // BRIDGE
          BRIDGE_sub(v, CircPool::BRIDGE_using_CX_0());
          done = true;
        }
      }
      if (!done) {  // default decomposition
        BRIDGE_sub(v, CircPool::BRIDGE_using_CX_1());
      }
    }
    return success;
  });
}

Transform decompose_CX_directed(const Architecture &arc) {
  return Transform([arc](Circuit &circ) {
    bool success = false;
    // Collect all CX type vertices
    std::vector<std::pair<Vertex, bool>> bin;
    for (Circuit::CommandIterator it = circ.begin(); it != circ.end(); ++it) {
      if (it->get_op_ptr()->get_type() == OpType::CX) {
        unit_vector_t qbs = it->get_args();
        node_vector_t nodes = {qbs.begin(), qbs.end()};
        if (!arc.edge_exists(nodes[0], nodes[1]) &&
            arc.edge_exists(nodes[1], nodes[0])) {
          // Implies CX gate is valid, and needs flipping to respect
          // Architecture
          bin.push_back({it.get_vertex(), false});
        }
      }
      if (it->get_op_ptr()->get_type() == OpType::Conditional) {
        const Conditional &b =
            static_cast<const Conditional &>(*it->get_op_ptr());
        if (b.get_op()->get_type() == OpType::CX) {
          qubit_vector_t qbs = it->get_qubits();
          node_vector_t nodes = {qbs.begin(), qbs.end()};
          if (!arc.edge_exists(nodes[0], nodes[1]) &&
              arc.edge_exists(nodes[1], nodes[0])) {
            // Implies CX gate is valid, and needs flipping to respect
            // Architecture
            bin.push_back({it.get_vertex(), true});
          }
        }
        if (b.get_op()->get_type() == OpType::CircBox) {
          qubit_vector_t qbs = it->get_qubits();
          std::shared_ptr<const Box> box_ptr =
              std::dynamic_pointer_cast<const Box>(b.get_op());
          qubit_vector_t all_qubits = box_ptr->to_circuit().get()->all_qubits();
          if (all_qubits.size() != 3)
            throw std::logic_error("Box being opened not a BRIDGE gate.");
          std::map<Qubit, Qubit> rmap = {
              {all_qubits[0], qbs[0]},
              {all_qubits[1], qbs[1]},
              {all_qubits[2], qbs[2]}};
          box_ptr->to_circuit()->rename_units(rmap);
          decompose_CX_directed(arc).apply(*box_ptr->to_circuit());
        }
      }
    }
    for (std::pair<Vertex, bool> v : bin) {
      if (!v.second) {
        Subcircuit sub = {
            circ.get_in_edges(v.first),
            circ.get_all_out_edges(v.first),
            {v.first}};
        circ.substitute(
            CircPool::CX_using_flipped_CX(), sub, Circuit::VertexDeletion::Yes);
      }
      if (v.second) {
        circ.substitute_conditional(
            CircPool::CX_using_flipped_CX(), v.first,
            Circuit::VertexDeletion::Yes);
      }
      success = true;
    }
    return success;
  });
}

Transform decompose_NPhasedX() {
  return Transform([](Circuit &circ) {
    VertexList bin;
    bool success = false;

    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::NPhasedX) {
        Gate_ptr g = as_gate_ptr(circ.get_Op_ptr_from_Vertex(v));
        unsigned n = g->n_qubits();
        Circuit sub(n);

        for (unsigned i = 0; i < n; ++i) {
          sub.add_op<unsigned>(OpType::PhasedX, g->get_params(), {i});
        }
        circ.substitute(sub, v, Circuit::VertexDeletion::No);
        bin.push_back(v);
        success = true;
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

///////////////////////////////
//     GlobalisePhasedX      //
///////////////////////////////

static unsigned n_distinct_beta(const Circuit &circ, const OptVertexVec &gates);
template <typename T>
static bool all_equal(const std::vector<T> &vs);

// Any PhasedX or NPhasedX gate is replaced by an NPhasedX that is global
Transform globalise_PhasedX(bool squash) {
  // The key bit: choose the decomposition strategy depending on the current
  // beta angles.
  //
  // Given the set of PhasedX gates to synthesise, choose which decomposition
  // strategy. Valid strategies are 0, 1, 2, corresponding to the number of
  // NPhasedX gates to be inserted.
  //
  // The current strategy is rather simple: it chooses to insert a single
  // NPhasedX whenever it would solve the current problem and there are further
  // PhasedX left (meaning that the rest of the computation can be deferred till
  // later), otherwise inserts 2x PhasedX.
  auto choose_strategy = [](const PhasedXFrontier &frontier, unsigned n_target,
                            unsigned n_all) {
    if (n_all == 0 || n_target == 0) {
      return 0;
    }
    if (n_target == 1 && n_all == 1) {
      return 1;
    }
    if (n_target == 1 && frontier.are_phasedx_left()) {
      return 1;
    }
    return 2;
  };

  // the actual transform
  return Transform([squash, choose_strategy](Circuit &circ) {
    // if we squash, we start by removing all NPhasedX gates
    if (squash) {
      Transforms::decompose_NPhasedX().apply(circ);
    }

    std::vector<unsigned> range_qbs(circ.n_qubits());
    std::iota(range_qbs.begin(), range_qbs.end(), 0);
    PhasedXFrontier frontier(circ);

    // find a total ordering of boundary gates (aka non-NPhasedX multi-qb gates)
    auto is_boundary = [&circ](Vertex v) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      return PhasedXFrontier::is_interval_boundary(op);
    };
    // get a lexicographic ordering of boundary vertices
    auto vertices_in_order = circ.vertices_in_order();
    auto r = vertices_in_order | boost::adaptors::filtered(is_boundary);
    OptVertexVec boundary_gates(r.begin(), r.end());
    // Add sentinel to process the gates after last boundary_gate.
    // This is required to force flush the frontier after the last multi-qubit
    // gate.
    OptVertex optv;
    boundary_gates.push_back(optv);

    // whether transform is successful (always true if squash=true)
    bool success = squash;

    // Loop through each multi-qb gate.
    // At each iteration, decide how to decompose the single-qb gates
    // immediately preceding the current multi-qb gate into global gates.
    //
    // The rationale: defer any computation until just before the next multi-qb
    // gate. At that point, the single-qb gates on the qubits where the multi-qb
    // gate is acting HAVE TO be decomposed. Figure out how to do that and take
    // care of the garbage you might have created on other qubits at a later
    // point.

    for (OptVertex v : boundary_gates) {
      // the qubits whose intervals must be decomposed into global gates
      std::set<unsigned> curr_qubits;
      while (true) {
        if (v) {
          curr_qubits = frontier.qubits_ending_in(*v);
        } else {
          curr_qubits.clear();
          for (unsigned i = 0; i < circ.n_qubits(); ++i) {
            curr_qubits.insert(curr_qubits.end(), i);
          }
        }
        if (squash) {
          frontier.squash_intervals();
        }
        OptVertexVec all_phasedx = frontier.get_all_beta_vertices();
        OptVertexVec curr_phasedx;
        for (unsigned q : curr_qubits) {
          curr_phasedx.push_back(all_phasedx[q]);
        }

        if (all_nullopt(curr_phasedx)) {
          // there is nothing to decompose anymore, move to next boundary_gate
          break;
        }
        if (all_equal(all_phasedx)) {
          // this is already a global NPhasedX gate, leave untouched
          frontier.skip_global_gates(1);
          continue;
        }

        // find best decomposition strategy
        unsigned n_curr_betas = n_distinct_beta(circ, curr_phasedx);
        unsigned n_all_betas = n_distinct_beta(circ, all_phasedx);
        unsigned strategy;
        if (squash) {
          strategy = choose_strategy(frontier, n_curr_betas, n_all_betas);
        } else {
          // if we don't squash we decompose each NPhasedX with 2x global
          strategy = 2;
        }
        switch (strategy) {
          case 0:
            // do nothing
            break;
          case 1:
            // insert one single global NPhasedX
            TKET_ASSERT(curr_qubits.size() > 0);
            frontier.insert_1_phasedx(*(curr_qubits.begin()));
            success = true;
            break;
          case 2:
            // insert two global NPhasedX
            frontier.insert_2_phasedx();
            success = true;
            break;
          default:
            throw NotValid("Invalid strategy in replace_non_global_phasedx");
        }
      }
      if (v) {
        frontier.next_multiqb(*v);
      }
    }
    TKET_ASSERT(frontier.is_finished());

    success |= absorb_Rz_NPhasedX().apply(circ);

    return success;
  });
}

// The number of distinct beta angles in `gates`.
//
// Performs O(n^2) pairwise equivalence comparisons between expressions to
// handle symbolic variables and floating point errors
unsigned n_distinct_beta(const Circuit &circ, const OptVertexVec &gates) {
  std::vector<Expr> vals;

  // collect all epxressions
  for (OptVertex v : gates) {
    if (v) {
      Expr angle = circ.get_Op_ptr_from_Vertex(*v)->get_params()[0];
      vals.push_back(angle);
    } else {
      vals.push_back(0.);
    }
  }

  unsigned n_distinct = 0;

  // perform pairwise equivalence checks
  for (unsigned i = 0; i < vals.size(); ++i) {
    bool is_unique = true;
    for (unsigned j = i + 1; j < vals.size(); ++j) {
      if (equiv_expr(vals[i], vals[j])) {
        is_unique = false;
        break;
      }
    }
    n_distinct += is_unique;
  }
  return n_distinct;
}

template <typename T>
bool all_equal(const std::vector<T> &vs) {
  if (vs.empty()) {
    return true;
  }
  T front = vs.front();
  for (auto v : vs) {
    if (front != v) {
      return false;
    }
  }
  return true;
}

}  // namespace Transforms

}  // namespace tket

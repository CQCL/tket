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

#include "tket/Transformations/Decomposition.hpp"

#include <functional>
#include <optional>
// replace with c++20 <ranges> when available
#include <boost/range/adaptor/filtered.hpp>
#include <stdexcept>

#include "tket/Architecture/Architecture.hpp"
#include "tket/Circuit/CircPool.hpp"
#include "tket/Circuit/CircUtils.hpp"
#include "tket/Converters/PhasePoly.hpp"
#include "tket/Gate/GatePtr.hpp"
#include "tket/OpType/OpType.hpp"
#include "tket/OpType/OpTypeFunctions.hpp"
#include "tket/OpType/OpTypeInfo.hpp"
#include "tket/Ops/OpPtr.hpp"
#include "tket/Transformations/BasicOptimisation.hpp"
#include "tket/Transformations/Rebase.hpp"
#include "tket/Transformations/Replacement.hpp"
#include "tket/Transformations/Transform.hpp"
#include "tket/Utils/Constants.hpp"
#include "tket/Utils/Expression.hpp"
#include "tket/Utils/MatrixAnalysis.hpp"

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
      Subcircuit sub = circ.singleton_subcircuit(v);
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
      Subcircuit sub = circ.singleton_subcircuit(v);
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
      Subcircuit sub = circ.singleton_subcircuit(v);
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
      Subcircuit sub = circ.singleton_subcircuit(v);
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
    bool conditional = false;
    while (optype == OpType::Conditional) {
      op = static_cast<const Conditional &>(*op).get_op();
      optype = op->get_type();
      conditional = true;
    }
    if (is_gate_type(optype) && !is_projective_type(optype) &&
        op->n_qubits() == 1 && optype != OpType::TK1) {
      std::vector<Expr> tk1_angs = as_gate_ptr(op)->get_tk1_angles();
      Circuit rep(1);
      rep.add_op<unsigned>(
          OpType::TK1, {tk1_angs[0], tk1_angs[1], tk1_angs[2]}, {0});
      if (conditional) {
        circ.substitute_conditional(rep, v, Circuit::VertexDeletion::No);
      } else {
        circ.substitute(rep, v, Circuit::VertexDeletion::No);
        circ.add_phase(tk1_angs[3]);
      }
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

static Circuit tk1_angles_to_circ(Expr a, Expr b, Expr c) {
  Circuit circ(1);

  std::optional<double> ea = eval_expr(a);
  std::optional<double> eb = eval_expr(b);
  std::optional<double> ec = eval_expr(c);

  // Remove global phases
  if (ea && *ea >= 2) {
    a -= 2;
    circ.add_phase(1);
  }
  if (eb && *eb >= 2) {
    b -= 2;
    circ.add_phase(1);
  }
  if (ec && *ec >= 2) {
    c -= 2;
    circ.add_phase(1);
  }

  if (!equiv_0(a, 2) || !equiv_0(b, 2) || !equiv_0(c, 2)) {
    circ.add_op<unsigned>(OpType::TK1, {c, b, a}, {0});
  }

  return circ;
}

Transform decompose_ZXZ_to_TK1() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    VertexVec inputs = circ.q_inputs();
    for (VertexVec::iterator i = inputs.begin(); i != inputs.end(); ++i) {
      Edge e = circ.get_nth_out_edge(*i, 0);
      Vertex v = circ.target(e);
      // Angles for current TK1
      std::array<Expr, 3> curr_angles;
      // Index from 0 <= curr_ind <= 2 into `curr_angles` collecting TK1 angles
      unsigned curr_ind = 0;
      // The beginning edge of a sequence of Rx and Rz gates
      std::optional<Edge> first_edge;
      // Vertices of the Rx/Rz sequence. Empty iff first_edge == std::nullopt
      VertexSet curr_vs;
      while (true) {
        Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
        bool repeat_loop = false;
        bool substitute_vertex = false;
        TKET_ASSERT(curr_ind < 3);
        switch (op->get_type()) {
          case OpType::Rz: {
            if (curr_ind == 1) {
              curr_angles[curr_ind++] = 0;
            }
            TKET_ASSERT(op->get_params().size() == 1);
            curr_angles[curr_ind++] = op->get_params().at(0);
            substitute_vertex = true;
            break;
          }
          case OpType::Rx: {
            if (curr_ind != 1) {
              curr_angles[curr_ind++] = 0;
            }
            TKET_ASSERT(op->get_params().size() == 1);
            if (curr_ind < 3) {
              curr_angles[curr_ind++] = op->get_params().at(0);
              substitute_vertex = true;
            } else {
              repeat_loop = true;
            }
            break;
          }
          default: {
            while (curr_ind < 3) {
              curr_angles[curr_ind++] = 0;
            }
          }
        }
        if (substitute_vertex) {
          curr_vs.insert(v);
          if (!first_edge) {
            first_edge = e;
          }
        }
        // Substitute sequence of RzRxRz with TK1
        auto substitute = [&]() {
          if (first_edge) {
            success = true;
            Circuit sub = tk1_angles_to_circ(
                curr_angles[0], curr_angles[1], curr_angles[2]);
            Subcircuit hole{{*first_edge}, {e}, {}, curr_vs};
            // Backup
            port_t backup_p = circ.get_target_port(e);
            // Substitute
            circ.substitute(sub, hole, Circuit::VertexDeletion::No);
            // Restore
            e = circ.get_nth_in_edge(v, backup_p);
            // Reset all tracking variables
            bin.insert(bin.end(), curr_vs.begin(), curr_vs.end());
            std::fill(curr_angles.begin(), curr_angles.end(), 0);
            curr_vs.clear();
            first_edge.reset();
          }
          curr_ind = 0;
        };
        // Depending on `substitute_vertex`, we either place the TK1 before
        // moving the edge forward or after
        if (curr_ind == 3 && !substitute_vertex) {
          substitute();
        }
        if (is_final_q_type(op->get_type())) {
          TKET_ASSERT(!substitute_vertex);
          break;
        }
        if (!repeat_loop) {
          // Move edge forward
          e = circ.get_next_edge(v, e);
          v = circ.target(e);
        }
        if (curr_ind == 3 && substitute_vertex) {
          substitute();
        }
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
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType optype = op->get_type();
      bool conditional = false;
      while (optype == OpType::Conditional) {
        op = static_cast<const Conditional &>(*op).get_op();
        optype = op->get_type();
        conditional = true;
      }
      if (optype == OpType::TK1) {
        success = true;
        bin.push_back(v);
        const std::vector<Expr> &params = op->get_params();
        Circuit newcirc =
            CircPool::tk1_to_rzrx(params[0], params[1], params[2]);
        if (conditional) {
          circ.substitute_conditional(newcirc, v, Circuit::VertexDeletion::No);
        } else {
          circ.substitute(newcirc, v, Circuit::VertexDeletion::No);
        }
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
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
        Subcircuit sub = {
            circ.get_in_edges(v), circ.get_linear_out_edges(v), {}, {}};
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
      std::vector<std::optional<Edge>> outs = circ.get_linear_out_edges(v);
      Vertex next = circ.target(*outs[0]);
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
              circ.get_nth_in_edge(last, 1) == *outs[1]) {
            // Recognise exp(-i XX * angle * pi/2)
            const Op_ptr op_ptr = get_op_ptr(OpType::XXPhase, angle);
            circ.dag[v] = {op_ptr};
            bin.push_back(next);
            circ.remove_vertex(
                next, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            bin.push_back(last);
            circ.remove_vertex(
                last, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            // Don't try to convert the other CX gate
            circ.set_vertex_Op_ptr(last, get_op_ptr(OpType::noop));
            circ.add_phase(phase);
            success = true;
            continue;
          }
        }
      }
      // Replace remaining CX gates
      Subcircuit sub = {circ.get_in_edges(v), outs, {}, {v}};
      bin.push_back(v);
      circ.substitute(
          CircPool::CX_using_XXPhase_1(), sub, Circuit::VertexDeletion::No);
      success = true;
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

/** Given TK2 angles, computes the fidelity that can be achieved using
 *  nb_cx CX gates.
 *
 *  @param k The TK2 angles.
 *  @param nb_cx The number of CX gates to be used for decomposition.
 *  @return The fidelity.
 */
static double get_CX_fidelity(const std::array<double, 3> &k, unsigned nb_cx) {
  TKET_ASSERT(nb_cx < 4);
  auto [a, b, c] = k;
  // gate fidelity achievable with 0,...,3 cnots
  // this is fully determined by the information content k and is optimal
  // see PhysRevA 71.062331 (2005) for more details on this
  switch (nb_cx) {
    case 0:
      return trace_fidelity(a, b, c);
    case 1:
      return trace_fidelity(0.5 - a, b, c);
    case 2:
      return trace_fidelity(0, 0, c);
    default:
      return 1.;
  }
}

// Try to decompose a TK2 gate using different gate sets, find the one with
// the highest fidelity.
// If no fidelities are provided, (best_optype, n_gates) is left unchanged.
static double best_noise_aware_decomposition(
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
        double gate_fid = std::visit(
            overloaded{
                // Constant value
                [](double arg) { return arg; },
                // A value depending on the angle
                [angles, n_zz](std::function<double(double)> arg) {
                  return (arg)(angles[n_zz - 1]);
                }},
            *fid.ZZPhase_fidelity);
        if (gate_fid < 0 || gate_fid > 1) {
          throw std::domain_error(
              "ZZPhase_fidelity returned a value outside of [0, 1].");
        }
        zz_fid *= gate_fid;
      }
      double nzz_fid = get_ZZPhase_fidelity(angles, n_zz) * zz_fid;
      // Use ZZPhase if fidelity is greater or it is equal but uses fewer gates
      if (nzz_fid - max_fid > EPS ||
          (nzz_fid - max_fid > -EPS && n_zz < n_gates)) {
        max_fid = nzz_fid;
        best_optype = OpType::ZZPhase;
        n_gates = n_zz;
      }
    }
  }

  return max_fid;
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
 * @param allow_swaps Whether implicit swaps are allowed.
 * @return Circuit TK2-equivalent circuit
 */
static Circuit TK2_replacement(
    std::array<Expr, 3> angles, const TwoQbFidelities &fid, bool allow_swaps) {
  if (!in_weyl_chamber(angles)) {
    throw std::domain_error("TK2 params are not normalised to Weyl chamber.");
  }
  OpType best_optype = OpType::CX;  // default to using CX
  unsigned n_gates = 3;             // default to 3x CX
  bool implicit_swap = false;       // default to no implicit swap

  // Only used when allow_swaps == true.
  Circuit pre, post;
  std::array<Expr, 3> angles_swapped;
  if (allow_swaps) {
    // Swapped circuit
    Circuit swap_circ(2);
    angles_swapped = angles;
    for (unsigned i = 0; i < 3; ++i) {
      angles_swapped[i] += 0.5;
    }
    std::tie(pre, angles_swapped, post) = normalise_TK2_angles(
        angles_swapped[0], angles_swapped[1], angles_swapped[2]);
    pre.add_phase(0.25);
  }

  // Try to evaluate exprs to doubles.
  std::array<double, 3> angles_eval;
  std::array<double, 3> angles_eval_swapped;
  unsigned last_angle = 0;
  for (; last_angle < 3; ++last_angle) {
    std::optional<double> eval = eval_expr_mod(angles[last_angle]);
    if (eval) {
      angles_eval[last_angle] = *eval;
    } else {
      break;
    }
    if (allow_swaps) {
      eval = eval_expr_mod(angles_swapped[last_angle]);
      TKET_ASSERT(eval);
      angles_eval_swapped[last_angle] = *eval;
    }
  }

  if (last_angle <= 2) {
    // Not all angles could be resolved numerically.
    // For symbolic angles, we can only provide an exact decomposition.
    best_exact_decomposition(angles, fid, best_optype, n_gates);
    if (allow_swaps) {
      OpType best_optype_swapped = OpType::CX;  // default to using CX
      unsigned n_gates_swapped = 3;             // default to 3x CX
      best_exact_decomposition(
          angles_swapped, fid, best_optype_swapped, n_gates_swapped);
      if (n_gates_swapped < n_gates) {
        n_gates = n_gates_swapped;
        best_optype = best_optype_swapped;
        angles = angles_swapped;
        implicit_swap = true;
      }
    }
  } else {
    // For non-symbolic angles, we can find the optimal number of gates
    // using the gate fidelities provided.
    double max_fid =
        best_noise_aware_decomposition(angles_eval, fid, best_optype, n_gates);
    if (allow_swaps) {
      OpType best_optype_swapped = OpType::CX;  // default to using CX
      unsigned n_gates_swapped = 3;             // default to 3x CX
      double max_fid_swapped = best_noise_aware_decomposition(
          angles_eval_swapped, fid, best_optype_swapped, n_gates_swapped);
      if (max_fid_swapped > max_fid ||
          ((max_fid_swapped - max_fid) < EPS && n_gates_swapped < n_gates)) {
        n_gates = n_gates_swapped;
        best_optype = best_optype_swapped;
        angles = angles_swapped;
        implicit_swap = true;
      }
    }
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
          TKET_ASSERT(!"Number of CX invalid in decompose_TK2");
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
          TKET_ASSERT(!"Number of ZZPhase invalid in decompose_TK2");
      }
      break;
    }
    default:
      throw BadOpType(
          "Unrecognised target OpType in decompose_TK2", best_optype);
  }

  if (implicit_swap) {
    Circuit swap(2);
    swap.add_op<unsigned>(OpType::SWAP, {0, 1});
    sub = pre >> sub >> post >> swap;
    sub.replace_SWAPs();
  }

  return sub;
}

Transform decompose_TK2(bool allow_swaps) {
  return decompose_TK2({}, allow_swaps);
}

Transform decompose_TK2(const TwoQbFidelities &fid, bool allow_swaps) {
  if (fid.ZZMax_fidelity) {
    if (*fid.ZZMax_fidelity < 0 || *fid.ZZMax_fidelity > 1) {
      throw std::domain_error("ZZMax fidelity must be between 0 and 1.");
    }
  }
  if (fid.CX_fidelity) {
    if (*fid.CX_fidelity < 0 || *fid.CX_fidelity > 1) {
      throw std::domain_error("CX fidelity must be between 0 and 1.");
    }
  }
  if (fid.ZZMax_fidelity && fid.ZZPhase_fidelity) {
    double ZZPhase_half = std::visit(
        overloaded{
            // A constant value.
            [](double arg) { return arg; },
            // A value depending on the input.
            [](std::function<double(double)> arg) { return (arg)(.5); }},
        *fid.ZZPhase_fidelity);
    if (*fid.ZZMax_fidelity < ZZPhase_half) {
      throw std::domain_error(
          "The ZZMax fidelity cannot be smaller than the ZZPhase(0.5) "
          "fidelity");
    }
  }
  return Transform([fid, allow_swaps](Circuit &circ) {
    bool success = false;

    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) != OpType::TK2) continue;

      success = true;
      auto params = circ.get_Op_ptr_from_Vertex(v)->get_params();
      TKET_ASSERT(params.size() == 3);
      std::array<Expr, 3> angles{params[0], params[1], params[2]};

      Circuit sub = TK2_replacement(angles, fid, allow_swaps);
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

Transform decompose_cliffords_std(bool tk2_to_cx) {
  return Transform([tk2_to_cx](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      OpType type = op->get_type();
      if (op->is_clifford()) {
        if (type != OpType::V && type != OpType::S && type != OpType::X &&
            type != OpType::Z && is_single_qubit_unitary_type(type)) {
          std::vector<Expr> tk1_param_exprs = as_gate_ptr(op)->get_tk1_angles();
          bool all_reduced = true;
          bool all_roundable = true;
          std::vector<int> iangles(3);
          for (int i = 0; i < 3; i++) {
            std::optional<double> reduced =
                eval_expr_mod(tk1_param_exprs[i], 4);
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
          Subcircuit sub = circ.singleton_subcircuit(v);
          bin.push_back(v);
          circ.substitute(replacement, sub, Circuit::VertexDeletion::No);
          circ.add_phase(tk1_param_exprs[3]);
          success = true;
        } else {
          switch (type) {
            case OpType::TK2: {
              if (tk2_to_cx) {
                auto params = op->get_params();
                TKET_ASSERT(params.size() == 3);
                Circuit replacement =
                    CircPool::TK2_using_CX(params[0], params[1], params[2]);
                decompose_cliffords_std().apply(replacement);
                bin.push_back(v);
                circ.substitute(replacement, v, Circuit::VertexDeletion::No);
                success = true;
              }
              break;
            }
            case OpType::NPhasedX: {
              auto params = op->get_params();
              unsigned n = circ.n_out_edges(v);
              Circuit replacement(n);
              for (unsigned i = 0; i < n; i++) {
                replacement.add_op<Qubit>(OpType::PhasedX, params, {Qubit(i)});
              }
              decompose_cliffords_std().apply(replacement);
              bin.push_back(v);
              circ.substitute(replacement, v, Circuit::VertexDeletion::No);
              success = true;
              break;
            }
            default:
              break;
          }
        }
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
                Subcircuit sub = circ.singleton_subcircuit(*it);
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

Transform decomp_boxes(
    const std::unordered_set<OpType> &excluded_types,
    const std::unordered_set<std::string> &excluded_opgroups,
    const std::optional<std::unordered_set<OpType>> &included_types,
    const std::optional<std::unordered_set<std::string>> &included_opgroups) {
  return Transform([=](Circuit &circ) {
    return circ.decompose_boxes_recursively(
        excluded_types, excluded_opgroups, included_types, included_opgroups);
  });
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
      std::vector<std::optional<Edge>> out_edges =
          circ.get_linear_out_edges(v.first);
      Subcircuit sub = {in_edges, out_edges, {}, {v.first}};

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
            {circ.get_target_port(*out_edges[0]),
             circ.get_target_port(*out_edges[1])});
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
        Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
        while (op->get_type() == OpType::Conditional) {
          op = static_cast<const Conditional &>(*op).get_op();
        }
        if (op->get_type() == OpType::BRIDGE) {
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
      std::vector<std::optional<Edge>> out_edges =
          circ.get_linear_out_edges(v.first);
      Subcircuit sub = {in_edges, out_edges, {}, {v.first}};

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
            circ.target(*out_edges[0]), circ.target(*out_edges[1]),
            circ.target(*out_edges[2])};
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
        Op_ptr op = it->get_op_ptr();
        while (op->get_type() == OpType::Conditional) {
          op = static_cast<const Conditional &>(*op).get_op();
        }
        if (op->get_type() == OpType::CX) {
          qubit_vector_t qbs = it->get_qubits();
          node_vector_t nodes = {qbs.begin(), qbs.end()};
          if (!arc.edge_exists(nodes[0], nodes[1]) &&
              arc.edge_exists(nodes[1], nodes[0])) {
            // Implies CX gate is valid, and needs flipping to respect
            // Architecture
            bin.push_back({it.get_vertex(), true});
          }
        }
        if (op->get_type() == OpType::CircBox) {
          qubit_vector_t qbs = it->get_qubits();
          std::shared_ptr<const Box> box_ptr =
              std::dynamic_pointer_cast<const Box>(op);
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
        Subcircuit sub = circ.singleton_subcircuit(v.first);
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

/* naive decomposition - there are cases we can do better if we can eg. ignore
 * phase */
Transform decomp_CCX() {
  return Transform([](Circuit &circ) {
    const Op_ptr ccx = get_op_ptr(OpType::CCX);
    return circ.substitute_all(CircPool::CCX_normal_decomp(), ccx);
  });
}

Transform decomp_controlled_Rys() {
  return Transform([](Circuit &circ) {
    bool success = decomp_CCX().apply(circ);
    auto [vit, vend] = boost::vertices(circ.dag);
    for (auto next = vit; vit != vend; vit = next) {
      ++next;
      Vertex v = *vit;
      const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
      unsigned arity = circ.n_in_edges(v);
      if (op->get_type() == OpType::CnRy) {
        success = true;
        Circuit rep = CircPool::CnRy_normal_decomp(op, arity);
        Subcircuit final_sub = circ.singleton_subcircuit(v);
        circ.substitute(rep, final_sub, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform decomp_arbitrary_controlled_gates() {
  static const std::set<OpType> cn_gate_set = {
      OpType::CCX, OpType::CnX,  OpType::CnRy, OpType::CnZ,
      OpType::CnY, OpType::CnRx, OpType::CnRz};
  std::set<OpType> all_gates;
  std::copy(
      all_gate_types().begin(), all_gate_types().end(),
      std::inserter(all_gates, all_gates.end()));
  OpTypeSet allowed_gate_set;
  std::set_difference(
      all_gates.begin(), all_gates.end(), cn_gate_set.begin(),
      cn_gate_set.end(),
      std::inserter(allowed_gate_set, allowed_gate_set.begin()));
  return rebase_factory(allowed_gate_set, CircPool::CX(), CircPool::tk1_to_tk1);
}

static void substitute_cnx(
    const Command &cmd, const std::optional<qubit_vector_t> &new_args,
    const bool &as_dagger, Circuit &circ,
    std::map<unsigned, Circuit> &cnx_cache) {
  Circuit cnx;
  qubit_vector_t original_args = cmd.get_qubits();
  bool conditional = cmd.get_op_ptr()->get_type() == OpType::Conditional;

  // Compute the decomposition
  auto cache = cnx_cache.find((unsigned)original_args.size());
  if (cache != cnx_cache.end()) {
    cnx = cache->second;
  } else {
    cnx = multi_controlled_to_2q(
        get_op_ptr(OpType::CnX, std::vector<Expr>(), original_args.size()));
    cnx_cache.insert({(unsigned)original_args.size(), cnx});
  }

  // Reorder the replacement circuit by renaming the qubits
  if (new_args != std::nullopt) {
    // The cnx circuit must only use default registers for this renaming to
    // work.
    TKET_ASSERT(cnx.is_simple());
    unit_map_t qmap;
    for (unsigned i = 0; i + 1 < new_args.value().size(); i++) {
      auto it =
          find(original_args.begin(), original_args.end(), new_args.value()[i]);
      int index = it - original_args.begin();
      qmap.insert({Qubit(i), Qubit(index)});
    }
    cnx.rename_units(qmap);
  }

  if (as_dagger) {
    cnx = cnx.dagger();
  }
  if (conditional) {
    circ.substitute_conditional(
        cnx, cmd.get_vertex(), Circuit::VertexDeletion::No);
  } else {
    circ.substitute(cnx, cmd.get_vertex(), Circuit::VertexDeletion::No);
  }
}

Transform cnx_pairwise_decomposition() {
  return Transform([](Circuit &circ) {
    bool success = commute_through_multis().apply(circ);
    success |= remove_redundancies().apply(circ);

    // Cache CnX decompositions
    std::map<unsigned, Circuit> cnx_cache;
    // Replaced vertices to delete at the end
    VertexList bin;

    // Find all CnX gates, including conditional ones
    std::vector<Command> commands = circ.get_commands();
    std::vector<int> cnx_indices;
    for (unsigned i = 0; i < commands.size(); i++) {
      Op_ptr op = commands[i].get_op_ptr();
      bool conditional = false;
      while (op->get_type() == OpType::Conditional) {
        op = static_cast<const Conditional &>(*op).get_op();
        conditional = true;
      }

      if (op->get_type() == OpType::CnX || op->get_type() == OpType::CCX) {
        if (op->n_qubits() < 3) {
          // Decompose as CX or X without further optimisation
          Circuit replacement = multi_controlled_to_2q(op);
          if (conditional) {
            circ.substitute_conditional(
                replacement, commands[i].get_vertex(),
                Circuit::VertexDeletion::No);
          } else {
            circ.substitute(
                replacement, commands[i].get_vertex(),
                Circuit::VertexDeletion::No);
          }
          bin.push_back(commands[i].get_vertex());
        } else {
          cnx_indices.push_back(i);
        }
      }
    }

    // If the number of CnX is odd, decompose the first one as it is
    unsigned itr = 0;
    if (cnx_indices.size() % 2 == 1) {
      substitute_cnx(
          commands[cnx_indices[0]], std::nullopt, false, circ, cnx_cache);
      bin.push_back(commands[cnx_indices[0]].get_vertex());
      itr = 1;
    }

    // For the rest of CnX, we reorder the control qubits of each pair
    // to maximise the chance of gate cancelling.
    for (; itr + 1 < cnx_indices.size(); itr += 2) {
      Command &cmd1 = commands[cnx_indices[itr]];
      Command &cmd2 = commands[cnx_indices[itr + 1]];

      qubit_vector_t args1 = cmd1.get_qubits();
      qubit_vector_t args2 = cmd2.get_qubits();

      qubit_vector_t new_args1 = args1;
      qubit_vector_t new_args2 = args2;

      bool cmd1_as_dagger = false;
      bool cmd2_as_dagger = false;

      if (args1.size() == 3 && args2.size() == 3) {
        // CCX gates are handled differently, there
        // is a maximum of 1 CX reduction (needs to be realised with
        // CliffordSimp), when they have args: [a,b,c] and [c,a,b]
        std::set<Qubit> args1_set(args1.begin(), args1.end());
        std::set<Qubit> args2_set(args2.begin(), args2.end());
        // Check if they are acting on the same set of qubits
        if (args1_set == args2_set) {
          // If they have the same target, we decompose them as CCX.dagger and
          // CCX
          if (new_args1.back() == new_args2.back()) {
            cmd1_as_dagger = true;
            new_args2 = new_args1;
          } else {
            // We find qubits a,b,c
            // such that the first CCX acts on [a,b,c]
            // and the second CCX acts on [c,a,b] without changing the semantics
            Qubit a = (args1[0] == args2.back()) ? args1[1] : args1[0];
            new_args1 = {a, args2.back(), args1.back()};
            new_args2 = {args1.back(), a, args2.back()};
          }
        }
        // If the two CCX gates are not acting on the same set of qubits
        // we don't change their arguments, and decompose them as they are.
      } else {
        // The more general case
        // we try to reorder the control args of the two gates to maximise
        // cancellation. The idea is to move the common control qubits to the
        // end of the control qubit lists e.g. C4X [q0, q1, q2, q3, q4] and C4X
        // [q1, q2, q4, q7, q6] will be reordered as C4X [q0, q3, q1, q2, q4]
        // and C4X [q4, q7, q1, q2, q6]
        qubit_vector_t common_ctrls;
        qubit_vector_t ctrl_args1 = args1;
        qubit_vector_t ctrl_args2 = args2;
        Qubit target1 = ctrl_args1.back();
        Qubit target2 = ctrl_args2.back();
        // Remove the target qubits so they only have control qubits
        ctrl_args1.pop_back();
        ctrl_args2.pop_back();
        // common_ctrls = ctrl_args1 intersect ctrl_args2
        std::sort(ctrl_args1.begin(), ctrl_args1.end());
        std::sort(ctrl_args2.begin(), ctrl_args2.end());
        std::set_intersection(
            ctrl_args1.begin(), ctrl_args1.end(), ctrl_args2.begin(),
            ctrl_args2.end(), std::back_inserter(common_ctrls));

        // reorder control qubits such that
        // new_args1 = (ctrl_args1 - common_ctrls) + common_ctrls + target1
        // new_args2 = (ctrl_args2 - common_ctrls) + common_ctrls + target2
        new_args1.clear();
        new_args2.clear();
        std::set_difference(
            ctrl_args1.begin(), ctrl_args1.end(), common_ctrls.begin(),
            common_ctrls.end(), std::back_inserter(new_args1));
        std::set_difference(
            ctrl_args2.begin(), ctrl_args2.end(), common_ctrls.begin(),
            common_ctrls.end(), std::back_inserter(new_args2));
        new_args1.insert(
            new_args1.end(), common_ctrls.begin(), common_ctrls.end());
        new_args2.insert(
            new_args2.end(), common_ctrls.begin(), common_ctrls.end());
        new_args1.push_back(target1);
        new_args2.push_back(target2);
        cmd2_as_dagger = true;
      }
      substitute_cnx(cmd1, new_args1, cmd1_as_dagger, circ, cnx_cache);
      substitute_cnx(cmd2, new_args2, cmd2_as_dagger, circ, cnx_cache);
      bin.push_back(cmd1.get_vertex());
      bin.push_back(cmd2.get_vertex());
    }

    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    success |= remove_redundancies().apply(circ);
    success |= !bin.empty();
    return success;
  });
}

}  // namespace Transforms

}  // namespace tket

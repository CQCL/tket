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
#include "Converters/PhasePoly.hpp"
#include "Gate/GatePtr.hpp"
#include "OpType/OpType.hpp"
#include "Ops/OpPtr.hpp"
#include "Replacement.hpp"
#include "Transform.hpp"
#include "Utils/Expression.hpp"

namespace tket {

static bool convert_to_zxz(Circuit &circ);
static bool convert_to_zyz(Circuit &circ);
static bool convert_to_xyx(Circuit &circ);

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
        !is_single_qubit_type(optype) && (optype != OpType::CX)) {
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
  bool success = (Transform::decompose_single_qubits_TK1() >>
                  Transform::decompose_tk1_to_rzrx())
                     .apply(circ);
  return success;
}

static bool convert_to_zyz(Circuit &circ) {
  static const Expr half = SymEngine::div(Expr(1), Expr(2));
  bool success = Transform::decompose_single_qubits_TK1().apply(circ);
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges(v) != 1) continue;
    const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    if (op->get_type() == OpType::tk1) {
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
  bool success = Transform::decompose_single_qubits_TK1().apply(circ);
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.n_in_edges(v) != 1) continue;
    const Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    if (op->get_type() == OpType::tk1) {
      std::vector<Expr> params = op->get_params();
      Circuit replacement(1);
      replacement.add_op<unsigned>(OpType::Ry, half, {0});
      replacement.add_op<unsigned>(OpType::Rx, params[2] + half, {0});
      replacement.add_op<unsigned>(OpType::Ry, params[1], {0});
      replacement.add_op<unsigned>(OpType::Rx, params[0] - half, {0});
      replacement.add_op<unsigned>(OpType::Ry, -half, {0});
      Transform::remove_redundancies().apply(replacement);
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

Transform Transform::decompose_multi_qubits_CX() {
  return Transform(convert_multiqs_CX);
}

static bool convert_singleqs_TK1(Circuit &circ) {
  bool success = false;
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    Op_ptr op = circ.get_Op_ptr_from_Vertex(v);
    OpType optype = op->get_type();
    if (is_single_qubit_type(optype) && !is_projective_type(optype) &&
        optype != OpType::tk1) {
      std::vector<Expr> tk1_angs = as_gate_ptr(op)->get_tk1_angles();
      Circuit rep(1);
      rep.add_op<unsigned>(
          OpType::tk1, {tk1_angs[0], tk1_angs[1], tk1_angs[2]}, {0});
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

Transform Transform::decompose_single_qubits_TK1() {
  return Transform(convert_singleqs_TK1);
}

Transform Transform::decompose_ZYZ_to_TK1() {
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
            circ.dag[v] = {get_op_ptr(OpType::tk1, new_params)};
          } else {
            circ.dag[v] = {get_op_ptr(OpType::tk1, {zero, zero, angle_1})};
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
          circ.dag[v] = {get_op_ptr(OpType::tk1, new_params)};
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

Transform Transform::decompose_ZX() { return Transform(convert_to_zxz); }

Transform Transform::decompose_ZY() { return Transform(convert_to_zyz); }

Transform Transform::decompose_XY() { return Transform(convert_to_xyx); }

Transform Transform::decompose_tk1_to_rzrx() {
  return Transform([](Circuit &circ) {
    bool success = false;
    auto [it, end] = boost::vertices(circ.dag);
    for (auto next = it; it != end; it = next) {
      ++next;
      if (circ.get_OpType_from_Vertex(*it) == OpType::tk1) {
        success = true;
        const Op_ptr g = circ.get_Op_ptr_from_Vertex(*it);
        const std::vector<Expr> &params = g->get_params();
        Circuit newcirc =
            Transform::tk1_to_rzrx(params[0], params[1], params[2]);
        Subcircuit sc = {
            {circ.get_in_edges(*it)}, {circ.get_all_out_edges(*it)}, {*it}};
        circ.substitute(newcirc, sc, Circuit::VertexDeletion::Yes);
      }
    }
    return success;
  });
}

Transform Transform::decompose_CX_to_ECR() {
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

Transform Transform::decompose_CX_to_HQS2() {
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
Transform Transform::decompose_ZX_to_HQS1() {
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
    Transform::remove_redundancies().apply(circ);
    return success;
  });
}

// Decompose CX into MolmerSorensen as:
// ---C---         -V-S-|-H-
//    |      -->    XX(pi/4)
// ---X---         -----|-Vdg-
Transform Transform::decompose_MolmerSorensen() {
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

Transform Transform::decompose_ZZPhase() {
  return Transform([](Circuit &circ) {
    bool success = Transform::decompose_PhaseGadgets().apply(circ);
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      if (circ.get_OpType_from_Vertex(v) == OpType::PhaseGadget) {
        const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
        circ.dag[v] = {get_op_ptr(OpType::ZZPhase, g->get_params()[0])};
      }
    }
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
 * The (i,j,k) entry in this table represents tk1(i/2, j/2, k/2).
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
 * Clifford circuit equivalent to tk1(i/2, j/2, k/2)
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

Transform Transform::decompose_cliffords_std() {
  return Transform([](Circuit &circ) {
    bool success = false;
    VertexList bin;
    BGL_FORALL_VERTICES(v, circ.dag, DAG) {
      OpType opT = circ.get_OpType_from_Vertex(v);
      if (opT == OpType::tk1 || opT == OpType::U3 || opT == OpType::U2 ||
          opT == OpType::U1 || opT == OpType::Rx || opT == OpType::Ry ||
          opT == OpType::Rz || opT == OpType::PhasedX) {
        const Op_ptr g = circ.get_Op_ptr_from_Vertex(v);
        std::vector<Expr> tk1_param_exprs = as_gate_ptr(g)->get_tk1_angles();
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
      }
    }
    circ.remove_vertices(
        bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
    return success;
  });
}

Transform Transform::decompose_ZX_to_cliffords() {
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

Transform Transform::decompose_PhaseGadgets() {
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
            (type == OpType::tk1 && equiv_0(g->get_params()[1]))) {
          Vertex last_v = circ.get_next_pair(next_v, outs[1]).first;
          if (circ.get_OpType_from_Vertex(last_v) == OpType::CX &&
              circ.get_nth_in_edge(last_v, 0) == outs[0]) {
            VertexList bin = {next_v, last_v};
            big_bin.push_back(next_v);
            big_bin.push_back(last_v);
            circ.remove_vertices(
                bin, Circuit::GraphRewiring::Yes, Circuit::VertexDeletion::No);
            Expr t = g->get_params()[0];
            if (type == OpType::tk1) {
              t += g->get_params()[2];
            }
            circ.dag[*it] = {get_op_ptr(OpType::PhaseGadget, {t}, 2)};
            if (type == OpType::U1) {
              circ.add_phase(t / 2);
            } else if (
                type == OpType::tk1 && equiv_val(g->get_params()[1], 2, 4)) {
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

Transform Transform::decomp_boxes() {
  return Transform([](Circuit &circ) { return circ.decompose_boxes(); });
}

Transform Transform::compose_phase_poly_boxes() {
  return Transform([](Circuit &circ) {
    CircToPhasePolyConversion conv = CircToPhasePolyConversion(circ);
    conv.convert();
    circ = conv.get_circuit();
    return true;
  });
}

Transform Transform::decompose_SWAP(const Circuit &replacement_circuit) {
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
  // swap_circ_1 leaves a CX{0,1} next to current CX{0,1}, if not we can assume
  // second case.
  if (port_comp == comp)
    circ.substitute(swap_circ_1, sub, Circuit::VertexDeletion::Yes);
  else
    circ.substitute(swap_circ_2, sub, Circuit::VertexDeletion::Yes);
}

Transform Transform::decompose_SWAP_to_CX(const Architecture &arc) {
  // Note that the default argument will be out of scope at call-time!
  //  => we replace the default empty Architecture with nullptr
  // we need to keep arc as a pointer as there is no such thing as
  // optional references in std
  const Architecture *arc_ptr = arc.n_uids() ? &arc : nullptr;
  return Transform([arc_ptr](Circuit &circ) {
    bool success = false;
    std::vector<std::pair<Vertex, bool>> bin;
    for (Circuit::CommandIterator it = circ.begin(); it != circ.end(); ++it) {
      if (it->get_op_ptr()->get_type() == OpType::SWAP) {
        unit_vector_t qbs = it->get_args();
        node_vector_t nodes = {qbs.begin(), qbs.end()};
        if (arc_ptr != nullptr && arc_ptr->uid_exists(nodes[0]) &&
            arc_ptr->uid_exists(nodes[1]) &&
            arc_ptr->connection_exists(nodes[1], nodes[0])) {
          bin.push_back({it.get_vertex(), true});
        } else {
          bin.push_back({it.get_vertex(), false});
        }
      }
    }

    for (std::pair<Vertex, bool> v : bin) {
      success = true;
      // Get predecessor vertices and successor vertices and find subcircuit for
      // replacement
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
        // savings. If SWAP doesn't lend itself to annihlation though, the SWAP
        // is inserted to reduce number of H gates added in a 'directed' CX
        // decomposition.
        // SWAP_using_CX_1 is added if the backwards direction is available on
        // the architecture
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

Transform Transform::decompose_BRIDGE_to_CX() {
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
      // Get predecessor vertices and successor vertices and find subcircuit for
      // replacement
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

Transform Transform::decompose_CX_directed(const Architecture &arc) {
  return Transform([arc](Circuit &circ) {
    bool success = false;
    // Collect all CX type vertices
    std::vector<std::pair<Vertex, bool>> bin;
    for (Circuit::CommandIterator it = circ.begin(); it != circ.end(); ++it) {
      if (it->get_op_ptr()->get_type() == OpType::CX) {
        unit_vector_t qbs = it->get_args();
        node_vector_t nodes = {qbs.begin(), qbs.end()};
        if (!arc.connection_exists(nodes[0], nodes[1]) &&
            arc.connection_exists(nodes[1], nodes[0])) {
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
          if (!arc.connection_exists(nodes[0], nodes[1]) &&
              arc.connection_exists(nodes[1], nodes[0])) {
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
          Transform::decompose_CX_directed(arc).apply(*box_ptr->to_circuit());
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
}  // namespace tket

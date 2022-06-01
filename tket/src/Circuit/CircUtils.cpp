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

#include "CircUtils.hpp"

#include <cmath>
#include <complex>
#include <vector>

#include "CircPool.hpp"
#include "Circuit/Circuit.hpp"
#include "Gate/GatePtr.hpp"
#include "Gate/GateUnitaryMatrixImplementations.hpp"
#include "Gate/Rotation.hpp"
#include "OpType/OpType.hpp"
#include "Ops/Op.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/Expression.hpp"
#include "Utils/MatrixAnalysis.hpp"
#include "Utils/UnitID.hpp"

namespace tket {

Eigen::Matrix2cd get_matrix(const Circuit &circ, const Vertex &vert) {
  const Op_ptr op = circ.get_Op_ptr_from_Vertex(vert);
  if (op->get_type() != OpType::TK1) {
    throw NotImplemented(
        "Cannot obtain matrix from gate: " + op->get_name() +
        ". Try rebasing to tket's internal representation.");
  }
  std::vector<Expr> ps = op->get_params();
  ps.push_back(0);
  return get_matrix_from_tk1_angles(ps);
}

Eigen::Matrix2cd get_matrix_from_circ(const Circuit &circ) {
  if (circ.n_qubits() != 1)
    throw CircuitInvalidity(
        "Getting Matrix: expected 1 qubit circuit, found " +
        std::to_string(circ.n_qubits()));
  Complex factor = std::exp(i_ * PI * eval_expr(circ.get_phase()).value());
  VertexVec qpath = circ.qubit_path_vertices(circ.all_qubits()[0]);
  unsigned N = qpath.size();
  if (N == 2) return factor * Eigen::Matrix2cd::Identity();
  Eigen::Matrix2cd m = get_matrix(circ, qpath[N - 2]);
  for (unsigned x = N - 3; x >= 1; --x) {
    m = m * get_matrix(circ, qpath[x]);
  }
  return factor * m;
}

Eigen::Matrix4cd get_matrix_from_2qb_circ(const Circuit &circ) {
  std::vector<QPathDetailed> all_paths = circ.all_qubit_paths();
  std::map<Vertex, Eigen::Matrix4cd> v_to_op;
  Eigen::Matrix4cd cnot, tonc, swap;
  // clang-format off
    cnot << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
    tonc << 1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0,
            0, 1, 0, 0;
    swap << 1, 0, 0, 0,
            0, 0, 1, 0,
            0, 1, 0, 0,
            0, 0, 0, 1;
  // clang-format on

  for (unsigned uqb = 0; uqb < 2; uqb++) {
    for (QPathDetailed::iterator it = all_paths[uqb].begin();
         it != all_paths[uqb].end(); ++it) {
      const Op_ptr o = circ.get_Op_ptr_from_Vertex(it->first);
      switch (o->get_type()) {
        case OpType::Input:
        case OpType::Create:
        case OpType::Output:
        case OpType::Discard: {
          v_to_op[it->first] = Eigen::Matrix4cd::Identity();
          break;
        }
        case OpType::SWAP: {
          if (uqb == 0) v_to_op[it->first] = swap;
          break;
        }
        case OpType::CX: {
          if (uqb == 0) {
            if (it->second == 0) {
              v_to_op[it->first] = cnot;
            } else {
              v_to_op[it->first] = tonc;
            }
          }
          break;
        }
        case OpType::TK2: {
          auto params = o->get_params();
          TKET_ASSERT(params.size() == 3);
          v_to_op[it->first] = get_matrix_from_2qb_circ(
              CircPool::TK2_using_CX(params[0], params[1], params[2]));
          break;
        }
        default: {
          if (o->get_desc().is_gate() && circ.n_in_edges(it->first) == 1 &&
              circ.n_out_edges(it->first) == 1) {
            const Op_ptr g = o;
            std::vector<Expr> ps = as_gate_ptr(g)->get_tk1_angles();
            Eigen::Matrix2cd mat = get_matrix_from_tk1_angles(ps);
            if (uqb == 0) {
              v_to_op[it->first] =
                  Eigen::kroneckerProduct(mat, Eigen::Matrix2cd::Identity());
            } else {
              v_to_op[it->first] =
                  Eigen::kroneckerProduct(Eigen::Matrix2cd::Identity(), mat);
            }
          } else
            throw NotImplemented(
                "Cannot obtain matrix from op " + o->get_name());
        }
      }
    }
  }
  Eigen::Matrix4cd m = Eigen::Matrix4cd::Identity();
  SliceVec slices = circ.get_slices();
  for (const Slice &s : slices) {
    for (const Vertex &v : s) {
      m = v_to_op[v] * m;
    }
  }
  return std::exp(i_ * PI * eval_expr(circ.get_phase()).value()) * m;
}

Circuit two_qubit_canonical(const Eigen::Matrix4cd &U) {
  if (!is_unitary(U)) {
    throw std::invalid_argument(
        "Non-unitary matrix passed to two_qubit_canonical");
  }

  auto [K1, A, K2] = get_information_content(U);

  K1 /= pow(K1.determinant(), 0.25);
  K2 /= pow(K2.determinant(), 0.25);
  auto [a, b, c] = A;

  // Decompose single qubits
  auto [K1a, K1b] = kronecker_decomposition(K1);
  auto [K2a, K2b] = kronecker_decomposition(K2);

  Circuit result(2);

  std::vector<double> angles_q0 = tk1_angles_from_unitary(K2a);
  std::vector<double> angles_q1 = tk1_angles_from_unitary(K2b);
  result.add_op<unsigned>(
      OpType::TK1, {angles_q0.begin(), angles_q0.end() - 1}, {0});
  result.add_op<unsigned>(
      OpType::TK1, {angles_q1.begin(), angles_q1.end() - 1}, {1});

  result.add_op<unsigned>(OpType::TK2, {a, b, c}, {0, 1});

  angles_q0 = tk1_angles_from_unitary(K1a);
  angles_q1 = tk1_angles_from_unitary(K1b);
  result.add_op<unsigned>(
      OpType::TK1, {angles_q0.begin(), angles_q0.end() - 1}, {0});
  result.add_op<unsigned>(
      OpType::TK1, {angles_q1.begin(), angles_q1.end() - 1}, {1});

  // this fixes phase if decomposition is exact
  Eigen::Matrix4cd reminder = get_matrix_from_2qb_circ(result).adjoint() * U;
  const Complex phase = reminder(0, 0);  // reminder = phase * I
  result.add_phase(arg(phase) / PI);
  return result;
}

// Factorize U as VD where V corresponds to a 2-CX circuit and
// D = diag(z, z*, z*, z). Return V and z.
static std::pair<Eigen::Matrix4cd, Complex> decompose_VD(
    const Eigen::Matrix4cd &U) {
  if (!is_unitary(U)) {
    throw std::invalid_argument("Non-unitary matrix passed to decompose_VD");
  }

  // The calculations below are derived from the proof of Proposition V.2 in
  // https://arxiv.org/abs/quant-ph/0308033.

  Eigen::Matrix4cd u = U / pow(U.determinant(), 0.25);
  Complex a = u(3, 0) * u(0, 3) - u(2, 0) * u(1, 3) - u(1, 0) * u(2, 3) +
              u(0, 0) * u(3, 3);
  Complex b = u(3, 1) * u(0, 2) - u(2, 1) * u(1, 2) - u(1, 1) * u(2, 2) +
              u(0, 1) * u(3, 2);
  // Now we want to find z such that |z|=1 and (az* - bz) is real.
  // The numerical stability of this function is a concern when a is close to
  // -b*. This problem can be demonstrated in artificially constructed
  // examples (passing unitaries very close to, but not quite, the identity to
  // the functions below). In these cases the product VD (or DV) may not
  // approximate U to within Eigen's default tolerance. Is there a way to
  // dodge this issue?
  Complex w = a + std::conj(b);
  double d = std::abs(w);
  // If w = 0 then we can set z = 1.
  Complex z = (d < EPS) ? 1. : w / d;
  Complex z0 = sqrt(z);
  Complex z1 = std::conj(z0);
  Eigen::Matrix4cd V = U;
  V.col(0) *= z1;
  V.col(1) *= z0;
  V.col(2) *= z0;
  V.col(3) *= z1;
  return {V, z0};
}

static void replace_TK2_2CX(Circuit &circ) {
  VertexList bin;
  BGL_FORALL_VERTICES(v, circ.dag, DAG) {
    if (circ.get_OpType_from_Vertex(v) != OpType::TK2) continue;
    auto params = circ.get_Op_ptr_from_Vertex(v)->get_params();
    TKET_ASSERT(params.size() == 3);
    TKET_ASSERT(equiv_0(params[2], 4));
    Circuit sub = CircPool::approx_TK2_using_2xCX(params[0], params[1]);
    bin.push_back(v);
    circ.substitute(sub, v, Circuit::VertexDeletion::No);
  }
  TKET_ASSERT(bin.size() == 1);
  circ.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);
}

std::pair<Circuit, Complex> decompose_2cx_VD(const Eigen::Matrix4cd &U) {
  auto [V, z0] = decompose_VD(U);
  Circuit circ = two_qubit_canonical(V);
  replace_TK2_2CX(circ);
  return {circ, z0};
}

std::pair<Circuit, Complex> decompose_2cx_DV(const Eigen::Matrix4cd &U) {
  auto [V, z0] = decompose_VD(U.adjoint());
  V.adjointInPlace();
  Circuit circ = two_qubit_canonical(V);
  replace_TK2_2CX(circ);
  return {circ, std::conj(z0)};
}

Circuit phase_gadget(unsigned n_qubits, const Expr &t, CXConfigType cx_config) {
  // Handle n_qubits==0 as a special case, or the calculations below
  // go badly wrong.
  Circuit new_circ(n_qubits);
  if (n_qubits == 0) {
    new_circ.add_phase(-t / 2);
    return new_circ;
  }
  switch (cx_config) {
    case CXConfigType::Snake: {
      for (unsigned i = n_qubits - 1; i != 0; --i) {
        unsigned j = i - 1;
        new_circ.add_op<unsigned>(OpType::CX, {i, j});
      }
      new_circ.add_op<unsigned>(OpType::Rz, t, {0});
      for (unsigned i = 0; i != n_qubits - 1; ++i) {
        unsigned j = i + 1;
        new_circ.add_op<unsigned>(OpType::CX, {j, i});
      }
      break;
    }
    case CXConfigType::Star: {
      for (unsigned i = n_qubits - 1; i != 0; --i) {
        new_circ.add_op<unsigned>(OpType::CX, {i, 0});
      }
      new_circ.add_op<unsigned>(OpType::Rz, t, {0});
      for (unsigned i = 1; i != n_qubits; ++i) {
        new_circ.add_op<unsigned>(OpType::CX, {i, 0});
      }
      break;
    }
    case CXConfigType::Tree: {
      unsigned complete_layers = floor(log2(n_qubits));
      unsigned dense_end = pow(2, complete_layers);
      for (unsigned i = 0; i < n_qubits - dense_end; i++)
        new_circ.add_op<unsigned>(
            OpType::CX, {dense_end + i, dense_end - 1 - i});
      for (unsigned step_size = 1; step_size < dense_end; step_size *= 2) {
        for (unsigned i = 0; i < dense_end; i += 2 * step_size)
          new_circ.add_op<unsigned>(OpType::CX, {i + step_size, i});
      }
      new_circ.add_op<unsigned>(OpType::Rz, t, {0});
      for (unsigned step_size = dense_end / 2; step_size >= 1; step_size /= 2) {
        for (unsigned i = 0; i < dense_end; i += 2 * step_size)
          new_circ.add_op<unsigned>(OpType::CX, {i + step_size, i});
      }
      for (unsigned i = 0; i < n_qubits - dense_end; i++)
        new_circ.add_op<unsigned>(
            OpType::CX, {dense_end + i, dense_end - 1 - i});
      break;
    }
    case CXConfigType::MultiQGate: {
      std::vector<std::vector<unsigned>> conjugations;
      int sign_correction = 1;
      for (int q = n_qubits - 1; q > 0; q -= 2) {
        if (q - 1 > 0) {
          unsigned i = q, j = q - 1;
          // this is only equal to the CX decompositions above
          // up to phase, but phase differences are cancelled out by
          // its dagger XXPhase(-1/2) below.
          new_circ.add_op<unsigned>(OpType::H, {i});
          new_circ.add_op<unsigned>(OpType::H, {j});
          new_circ.add_op<unsigned>(OpType::XXPhase3, 0.5, {i, j, 0});
          sign_correction *= -1;
          conjugations.push_back({i, j, 0});
        } else {
          unsigned i = q;
          new_circ.add_op<unsigned>(OpType::CX, {i, 0});
          conjugations.push_back({i, 0});
        }
      }
      new_circ.add_op<unsigned>(OpType::Rz, sign_correction * t, {0});
      for (const auto &conj : conjugations) {
        if (conj.size() == 2) {
          new_circ.add_op<unsigned>(OpType::CX, conj);
        } else {
          TKET_ASSERT(conj.size() == 3);
          new_circ.add_op<unsigned>(OpType::XXPhase3, -0.5, conj);
          new_circ.add_op<unsigned>(OpType::H, {conj[0]});
          new_circ.add_op<unsigned>(OpType::H, {conj[1]});
        }
      }
      break;
    }
  }
  return new_circ;
}

Circuit pauli_gadget(
    const std::vector<Pauli> &paulis, const Expr &t, CXConfigType cx_config) {
  unsigned n = paulis.size();
  Circuit circ(n);
  std::vector<unsigned> qubits;
  for (unsigned i = 0; i < n; i++) {
    switch (paulis[i]) {
      case Pauli::I:
        break;
      case Pauli::X:
        circ.add_op<unsigned>(OpType::H, {i});
        qubits.push_back(i);
        break;
      case Pauli::Y:
        circ.add_op<unsigned>(OpType::V, {i});
        qubits.push_back(i);
        break;
      case Pauli::Z:
        qubits.push_back(i);
        break;
    }
  }
  if (qubits.empty()) {
    circ.add_phase(-t / 2);
  } else {
    Vertex v = circ.add_op<unsigned>(OpType::PhaseGadget, t, qubits);
    Circuit cx_gadget = phase_gadget(circ.n_in_edges(v), t, cx_config);
    Subcircuit sub = {circ.get_in_edges(v), circ.get_all_out_edges(v), {v}};
    circ.substitute(cx_gadget, sub, Circuit::VertexDeletion::Yes);
    for (unsigned i = 0; i < n; i++) {
      switch (paulis[i]) {
        case Pauli::I:
          break;
        case Pauli::X:
          circ.add_op<unsigned>(OpType::H, {i});
          break;
        case Pauli::Y:
          circ.add_op<unsigned>(OpType::Vdg, {i});
          break;
        case Pauli::Z:
          break;
      }
    }
  }
  return circ;
}

void replace_CX_with_TK2(Circuit &c) {
  static const Op_ptr cx = std::make_shared<Gate>(OpType::CX);
  c.substitute_all(CircPool::CX_using_TK2(), cx);
}

Circuit with_TK2(Gate_ptr op) {
  std::vector<Expr> params = op->get_params();
  unsigned n = op->n_qubits();
  if (n == 0) {
    return Circuit();
  } else if (n == 1) {
    Circuit c(1);
    c.add_op(op, std::vector<unsigned>{0});
    return c;
  } else if (n == 2 && op->free_symbols().empty()) {
    Eigen::Matrix4cd U = op->get_unitary();
    auto [K1, A, K2] = get_information_content(U);
    // Decompose single qubits
    auto [K1a, K1b] = kronecker_decomposition(K1);
    auto [K2a, K2b] = kronecker_decomposition(K2);
    Circuit c(2);
    std::vector<double> angles_K1a = tk1_angles_from_unitary(K1a);
    std::vector<double> angles_K1b = tk1_angles_from_unitary(K1b);
    std::vector<double> angles_K2a = tk1_angles_from_unitary(K2a);
    std::vector<double> angles_K2b = tk1_angles_from_unitary(K2b);
    c.add_op<unsigned>(
        OpType::TK1, {angles_K2a.begin(), angles_K2a.end() - 1}, {0});
    c.add_op<unsigned>(
        OpType::TK1, {angles_K2b.begin(), angles_K2b.end() - 1}, {1});
    double alpha = std::get<0>(A);
    double beta = std::get<1>(A);
    double gamma = std::get<2>(A);
    c.add_op<unsigned>(OpType::TK2, {alpha, beta, gamma}, {0, 1});
    c.add_op<unsigned>(
        OpType::TK1, {angles_K1a.begin(), angles_K1a.end() - 1}, {0});
    c.add_op<unsigned>(
        OpType::TK1, {angles_K1b.begin(), angles_K1b.end() - 1}, {1});

    // Correct phase by computing the unitary and comparing with U:
    Eigen::Matrix4cd V_K1 = Eigen::KroneckerProduct(
        get_matrix_from_tk1_angles(
            {angles_K1a[0], angles_K1a[1], angles_K1a[2], 0}),
        get_matrix_from_tk1_angles(
            {angles_K1b[0], angles_K1b[1], angles_K1b[2], 0}));
    Eigen::Matrix4cd V_A =
        internal::GateUnitaryMatrixImplementations::TK2(alpha, beta, gamma);
    Eigen::Matrix4cd V_K2 = Eigen::KroneckerProduct(
        get_matrix_from_tk1_angles(
            {angles_K2a[0], angles_K2a[1], angles_K2a[2], 0}),
        get_matrix_from_tk1_angles(
            {angles_K2b[0], angles_K2b[1], angles_K2b[2], 0}));
    Eigen::Matrix4cd V = V_K1 * V_A * V_K2;
    Eigen::Matrix4cd R = V.adjoint() * U;
    const Complex phase = R(0, 0);  // R = phase * I
    c.add_phase(arg(phase) / PI);

    return c;
  }
  // Now the non-trivial cases.
  switch (op->get_type()) {
    case OpType::ISWAP:
      return CircPool::ISWAP_using_TK2(params[0]);
    case OpType::PhasedISWAP:
      return CircPool::PhasedISWAP_using_TK2(params[0], params[1]);
    case OpType::XXPhase:
      return CircPool::XXPhase_using_TK2(params[0]);
    case OpType::YYPhase:
      return CircPool::YYPhase_using_TK2(params[0]);
    case OpType::ZZPhase:
      return CircPool::ZZPhase_using_TK2(params[0]);
    case OpType::NPhasedX:
      return CircPool::NPhasedX_using_PhasedX(n, params[0], params[1]);
    case OpType::ESWAP:
      return CircPool::ESWAP_using_TK2(params[0]);
    case OpType::FSim:
      return CircPool::FSim_using_TK2(params[0], params[1]);
    case OpType::CRx:
      return CircPool::CRx_using_TK2(params[0]);
    case OpType::CRy:
      return CircPool::CRy_using_TK2(params[0]);
    case OpType::CRz:
      return CircPool::CRz_using_TK2(params[0]);
    case OpType::CU1:
      return CircPool::CU1_using_TK2(params[0]);
    case OpType::XXPhase3:
      return CircPool::XXPhase3_using_TK2(params[0]);
    case OpType::CCX:
    case OpType::CSWAP:
    case OpType::BRIDGE:
    case OpType::CU3:
    case OpType::PhaseGadget: {
      // As a first, inefficient, solution, decompose these into CX and then
      // replace each CX with a TK2 (and some single-qubit gates).
      // TODO Find more efficient decompositions for these gates.
      Circuit c = with_CX(op);
      replace_CX_with_TK2(c);
      return c;
    }
    default:
      throw CircuitInvalidity("Cannot decompose " + op->get_name());
  }
}

Circuit with_CX(Gate_ptr op) {
  OpType optype = op->get_type();
  std::vector<Expr> params = op->get_params();
  unsigned n = op->n_qubits();
  if (n == 0) {
    return Circuit();
  } else if (n == 1) {
    Circuit c(1);
    c.add_op(op, std::vector<unsigned>{0});
    return c;
  }
  switch (optype) {
    case OpType::CX: {
      Circuit c(2);
      c.add_op(op, std::vector<unsigned>{0, 1});
      return c;
    }
    case OpType::CCX:
      return CircPool::CCX_normal_decomp();
    case OpType::CY:
      return CircPool::CY_using_CX();
    case OpType::CZ:
      return CircPool::CZ_using_CX();
    case OpType::CH:
      return CircPool::CH_using_CX();
    case OpType::CV:
      return CircPool::CV_using_CX();
    case OpType::CVdg:
      return CircPool::CVdg_using_CX();
    case OpType::CSX:
      return CircPool::CSX_using_CX();
    case OpType::CSXdg:
      return CircPool::CSXdg_using_CX();
    case OpType::CRz:
      return CircPool::CRz_using_CX(params[0]);
    case OpType::CRx:
      return CircPool::CRx_using_CX(params[0]);
    case OpType::CRy:
      return CircPool::CRy_using_CX(params[0]);
    case OpType::CU1:
      return CircPool::CU1_using_CX(params[0]);
    case OpType::CU3:
      return CircPool::CU3_using_CX(params[0], params[1], params[2]);
    case OpType::PhaseGadget:
      return phase_gadget(n, params[0], CXConfigType::Snake);
    case OpType::SWAP:
      return CircPool::SWAP_using_CX_0();
    case OpType::CSWAP:
      return CircPool::CSWAP_using_CX();
    case OpType::BRIDGE:
      return CircPool::BRIDGE_using_CX_0();
    case OpType::ECR:
      return CircPool::ECR_using_CX();
    case OpType::ISWAP:
      return CircPool::ISWAP_using_CX(params[0]);
    case OpType::ZZMax:
      return CircPool::ZZMax_using_CX();
    case OpType::XXPhase:
      return CircPool::XXPhase_using_CX(params[0]);
    case OpType::YYPhase:
      return CircPool::YYPhase_using_CX(params[0]);
    case OpType::ZZPhase:
      return CircPool::ZZPhase_using_CX(params[0]);
    case OpType::TK2:
      return CircPool::TK2_using_CX(params[0], params[1], params[2]);
    case OpType::XXPhase3:
      return CircPool::XXPhase3_using_CX(params[0]);
    case OpType::ESWAP:
      return CircPool::ESWAP_using_CX(params[0]);
    case OpType::FSim:
      return CircPool::FSim_using_CX(params[0], params[1]);
    case OpType::Sycamore:
      return CircPool::FSim_using_CX(1. / 2., 1. / 6.);
    case OpType::ISWAPMax:
      return CircPool::ISWAP_using_CX(1.);
    case OpType::PhasedISWAP:
      return CircPool::PhasedISWAP_using_CX(params[0], params[1]);
    case OpType::NPhasedX:
      return CircPool::NPhasedX_using_PhasedX(n, params[0], params[1]);
    default:
      throw CircuitInvalidity("Cannot decompose " + op->get_name());
  }
}

#define CNXTYPE(n) \
  (((n) == 2) ? OpType::CX : ((n) == 3) ? OpType::CCX : OpType::CnX)

/**
 * Construct a circuit representing CnU1 using U1, CX, CCX and CnX gates.
 */
static Circuit CnU1(unsigned n_controls, Expr lambda) {
  // CnU1(x) decomposes recursively as:
  // C{n-1}U1(x/2)[ctrls]; U1(x/2)[tgt]; CnX; U1(-x/2)[tgt]; CnX
  // We don't actually use recursion; just iterate starting with the first U1:
  Circuit c(n_controls + 1);
  Expr x = lambda / (1u << n_controls);
  c.add_op<unsigned>(OpType::U1, x, {0});
  std::vector<unsigned> cnx_qbs = {0};
  for (unsigned i = 0; i < n_controls; i++) {
    cnx_qbs.push_back(i + 1);
    c.add_op<unsigned>(OpType::U1, x, {i + 1});
    c.add_op<unsigned>(CNXTYPE(i + 2), cnx_qbs);
    c.add_op<unsigned>(OpType::U1, -x, {i + 1});
    c.add_op<unsigned>(CNXTYPE(i + 2), cnx_qbs);
    x *= 2;
  }
  return c;
}

Circuit with_controls(const Circuit &c, unsigned n_controls) {
  if (c.n_bits() != 0 || !c.is_simple()) {
    throw CircuitInvalidity("Only default qubit register allowed");
  }
  if (c.has_implicit_wireswaps()) {
    throw CircuitInvalidity("Circuit has implicit wireswaps");
  }

  // Dispose of the trivial case
  if (n_controls == 0) {
    return c;
  }

  unsigned c_n_qubits = c.n_qubits();

  // 1. Rebase to {CX, CCX, CnX, CnRy} and single-qubit gates
  Circuit c1(c);
  VertexList bin;
  BGL_FORALL_VERTICES(v, c1.dag, DAG) {
    Op_ptr op = c1.get_Op_ptr_from_Vertex(v);
    OpType optype = op->get_type();
    if (is_gate_type(optype)) {
      if (is_projective_type(optype)) {
        throw CircuitInvalidity("Projective operations present");
      }
      if (is_box_type(optype)) {
        throw CircuitInvalidity("Undecomposed boxes present");
      }
      if (is_single_qubit_type(optype)) {
        continue;
      }
      if (optype == OpType::CX || optype == OpType::CCX ||
          optype == OpType::CnX || optype == OpType::CnRy) {
        continue;
      }
      Circuit replacement = with_CX(as_gate_ptr(op));
      c1.substitute(replacement, v, Circuit::VertexDeletion::No);
      bin.push_back(v);
    }
  }
  c1.remove_vertices(
      bin, Circuit::GraphRewiring::No, Circuit::VertexDeletion::Yes);

  // Capture the phase. We may adjust this during replacements below.
  Expr a = c1.get_phase();

  // 2. Replace all gates with controlled versions
  Circuit c2(n_controls + c_n_qubits);
  for (Circuit::CommandIterator cit = c1.begin(); cit != c1.end(); ++cit) {
    Op_ptr op = cit->get_op_ptr();
    unit_vector_t args = cit->get_args();
    unsigned n_args = args.size();
    unsigned n_new_args = n_controls + n_args;
    qubit_vector_t new_args(n_new_args);
    for (unsigned i = 0; i < n_controls; i++) {
      new_args[i] = Qubit(i);
    }
    for (unsigned i = 0; i < n_args; i++) {
      new_args[n_controls + i] = Qubit(n_controls + args[i].index()[0]);
    }
    OpType optype = op->get_type();
    std::vector<Expr> params = op->get_params();
    switch (optype) {
      case OpType::noop:
        break;
      case OpType::X:
      case OpType::CX:
      case OpType::CCX:
      case OpType::CnX:
        c2.add_op<Qubit>(CNXTYPE(n_new_args), new_args);
        break;
      case OpType::Ry:
      case OpType::CnRy:
        c2.add_op<Qubit>(OpType::CnRy, params, new_args);
        break;
      default: {
        std::vector<Expr> tk1_angles = as_gate_ptr(op)->get_tk1_angles();
        Expr theta = tk1_angles[1];
        Expr phi = tk1_angles[0] - 0.5;
        Expr lambda = tk1_angles[2] + 0.5;
        Expr t = tk1_angles[3] - 0.5 * (tk1_angles[0] + tk1_angles[2]);
        // Operation is U3(theta, phi, lambda) + phase t.
        // First absorb t in the overall phase.
        a += t;
        // Construct a multi-controlled U3, by extending the standard
        // CU3-to-CX decomposition.
        Qubit target = new_args[n_controls];
        Circuit cnu1 = CnU1(n_controls - 1, 0.5 * (lambda + phi));
        c2.append(cnu1);
        c2.add_op<Qubit>(OpType::U1, 0.5 * (lambda - phi), {target});
        c2.add_op<Qubit>(CNXTYPE(n_new_args), new_args);
        c2.add_op<Qubit>(
            OpType::U3, {-0.5 * theta, 0, -0.5 * (lambda + phi)}, {target});
        c2.add_op<Qubit>(CNXTYPE(n_new_args), new_args);
        c2.add_op<Qubit>(OpType::U3, {0.5 * theta, phi, 0}, {target});
      } break;
    }
  }

  // 3. Account for phase by appending a CnU1 to the control qubits.
  if (!equiv_0(a)) {
    Circuit cnu1 = CnU1(n_controls - 1, a);
    c2.append(cnu1);
  }

  return c2;
}

#undef CNXTYPE

}  // namespace tket

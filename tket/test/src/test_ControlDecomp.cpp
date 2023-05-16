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

#include <boost/dynamic_bitset.hpp>
#include <catch2/catch_test_macros.hpp>
#include <numeric>

#include "Simulation/ComparisonFunctions.hpp"
#include "testutil.hpp"
#include "tket/Circuit/CircPool.hpp"
#include "tket/Gate/GateUnitaryMatrix.hpp"
#include "tket/Gate/SymTable.hpp"
#include "tket/Simulation/CircuitSimulator.hpp"
#include "tket/Transformations/CliffordReductionPass.hpp"
#include "tket/Transformations/Decomposition.hpp"
#include "tket/Transformations/OptimisationPass.hpp"
#include "tket/Transformations/Replacement.hpp"
#include "tket/Transformations/Transform.hpp"

namespace tket {
namespace test_ControlDecomp {

static bool approx_equal(const Complex& c1, const Complex& c2) {
  return (std::abs(c1 - c2) < ERR_EPS);
}

static bool check_incrementer_borrow_n_qubits(const unsigned n) {
  Circuit inc = CircPool::incrementer_borrow_n_qubits(n);
  bool correct = true;
  const StateVector sv = tket_sim::get_statevector(inc);
  for (unsigned i = 0; i < sv.size(); ++i) {
    // incremented the |0...00> state to be |0...10> incl. garbage qubits
    // (depending on def. of qubit significance)
    if (i == pow(2, 2 * n - 2))
      correct &= (std::abs(sv[i]) > EPS);
    else
      correct &= (std::abs(sv[i]) < ERR_EPS);
  }

  Circuit xcirc = Circuit(2 * n);
  for (unsigned i = 1; i < 2 * n; i += 2)
    xcirc.add_op<unsigned>(OpType::X, {i});
  xcirc.append(inc);
  const StateVector sv2 = tket_sim::get_statevector(xcirc);
  for (unsigned i = 0; i < sv2.size(); ++i) {
    if (i == 0)
      correct &= (std::abs(sv2[i]) > EPS);
    else
      correct &= (std::abs(sv2[i]) < ERR_EPS);
  }
  return correct;
}

static bool check_incrementer_borrow_1_qubit(const unsigned n) {
  Circuit inc = CircPool::incrementer_borrow_1_qubit(n);
  REQUIRE(inc.n_vertices() - inc.n_gates() == (n + 1) * 2);
  Transforms::synthesise_tket().apply(inc);
  const StateVector sv = tket_sim::get_statevector(inc);
  bool correct = true;
  for (unsigned i = 0; i < sv.size(); ++i) {
    // |00...0> -> |00...1>
    if (i == pow(2, n))
      correct &= (std::abs(sv[i]) > EPS);
    else
      correct &= (std::abs(sv[i]) < ERR_EPS);
  }
  Circuit xcirc = Circuit(n + 1);
  for (unsigned i = 0; i < n; i++) xcirc.add_op<unsigned>(OpType::X, {i});
  xcirc.append(inc);
  const StateVector sv2 = tket_sim::get_statevector(xcirc);
  for (unsigned i = 0; i < sv2.size(); ++i) {
    // |01...1> -> |00...0>
    if (i == 0)
      correct &= (std::abs(sv2[i]) > EPS);
    else
      correct &= (std::abs(sv2[i]) < ERR_EPS);
  }
  return correct;
}

static bool check_incrementer_linear_depth(
    const unsigned n, const unsigned number) {
  boost::dynamic_bitset<> in_bits(n, number);
  Circuit circ(n);
  for (unsigned i = 0; i < n; i++) {
    if (in_bits[i]) {
      circ.add_op<unsigned>(OpType::X, {i});
    }
  }
  Circuit inc = CircPool::incrementer_linear_depth(n);
  circ.append(inc);

  unsigned long correct_out_long = in_bits.to_ulong() + 1UL;
  boost::dynamic_bitset<> correct_out_bits(n, correct_out_long);

  // Get the index of the entry that should have a "1" (with a phase difference)
  unsigned sv_set_idx = 0;
  for (unsigned i = 0; i < correct_out_bits.size(); i++) {
    if (correct_out_bits[i] == 1) {
      sv_set_idx = sv_set_idx + (1 << (n - i - 1));
    }
  }
  const StateVector sv = tket_sim::get_statevector(circ);
  for (unsigned i = 0; i < sv.size(); ++i) {
    if (i == sv_set_idx) {
      if (std::abs(std::abs(sv[i]) - 1) >= ERR_EPS) return false;
    } else {
      if (std::abs(sv[i]) >= ERR_EPS) return false;
    }
  }
  return true;
}

// Explicitly construct CnU matrix
static Eigen::MatrixXcd get_CnU_matrix(
    unsigned n_controls, const Eigen::Matrix2cd& U) {
  unsigned m_size = pow(2, n_controls + 1);
  Eigen::MatrixXcd correct_matrix = Eigen::MatrixXcd::Identity(m_size, m_size);
  correct_matrix(m_size - 2, m_size - 2) = U(0, 0);
  correct_matrix(m_size - 2, m_size - 1) = U(0, 1);
  correct_matrix(m_size - 1, m_size - 2) = U(1, 0);
  correct_matrix(m_size - 1, m_size - 1) = U(1, 1);
  return correct_matrix;
}

static Eigen::MatrixXcd get_CnX_matrix(unsigned n_controls) {
  Eigen::Matrix2cd x = GateUnitaryMatrix::get_unitary(OpType::X, 1, {});
  return get_CnU_matrix(n_controls, x);
}
static Eigen::MatrixXcd get_CnY_matrix(unsigned n_controls) {
  Eigen::Matrix2cd y = GateUnitaryMatrix::get_unitary(OpType::Y, 1, {});
  return get_CnU_matrix(n_controls, y);
}
static Eigen::MatrixXcd get_CnZ_matrix(unsigned n_controls) {
  Eigen::Matrix2cd z = GateUnitaryMatrix::get_unitary(OpType::Z, 1, {});
  return get_CnU_matrix(n_controls, z);
}

SCENARIO("Test decomposition using CX") {
  OpType cntype;
  std::function<Eigen::MatrixXcd(unsigned)> matrix_func;
  WHEN("CnX") {
    cntype = OpType::CnX;
    matrix_func = get_CnX_matrix;
  }
  WHEN("CnY") {
    cntype = OpType::CnY;
    matrix_func = get_CnY_matrix;
  }
  WHEN("CnZ") {
    cntype = OpType::CnZ;
    matrix_func = get_CnZ_matrix;
  }
  std::vector<std::pair<unsigned, unsigned>> n_ctr_2q_count{
      {3, 14}, {4, 36}, {6, 120}};
  for (auto pair : n_ctr_2q_count) {
    const Op_ptr op = get_op_ptr(cntype, std::vector<Expr>(), pair.first + 1);
    Circuit decomposed_circ = CX_circ_from_multiq(op);
    auto u = tket_sim::get_unitary(decomposed_circ);
    REQUIRE((matrix_func(pair.first) - u).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(decomposed_circ.count_gates(OpType::CX) == pair.second);
  }
}

SCENARIO("Test decomposition using TK2") {
  OpType cntype;
  std::function<Eigen::MatrixXcd(unsigned)> matrix_func;
  WHEN("CnX") {
    cntype = OpType::CnX;
    matrix_func = get_CnX_matrix;
  }
  WHEN("CnY") {
    cntype = OpType::CnY;
    matrix_func = get_CnY_matrix;
  }
  WHEN("CnZ") {
    cntype = OpType::CnZ;
    matrix_func = get_CnZ_matrix;
  }
  std::vector<std::pair<unsigned, unsigned>> n_ctr_2q_count{
      {3, 14}, {4, 36}, {6, 61}};
  for (auto pair : n_ctr_2q_count) {
    const Op_ptr op = get_op_ptr(cntype, std::vector<Expr>(), pair.first + 1);
    Circuit decomposed_circ = TK2_circ_from_multiq(op);
    auto u = tket_sim::get_unitary(decomposed_circ);
    REQUIRE((matrix_func(pair.first) - u).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(decomposed_circ.count_gates(OpType::TK2) == pair.second);
  }
}

SCENARIO("Decompose some circuits with CCX gates") {
  GIVEN("Two CCX gates") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    Circuit circ2(3);
    const StateVector sv2 = tket_sim::get_statevector(circ2);
    Transforms::decomp_CCX().apply(circ);
    const StateVector sv1 = tket_sim::get_statevector(circ);
    REQUIRE(tket_sim::compare_statevectors_or_unitaries(sv1, sv2));

    WHEN("Check gate numbering") {
      Circuit circ3(3);
      circ3.add_op<unsigned>(OpType::CCX, {0, 1, 2});
      Transforms::decomp_CCX().apply(circ3);
      REQUIRE(circ3.n_gates() == 15);
      REQUIRE(circ3.n_vertices() == 21);
      REQUIRE(circ3.n_qubits() == 3);
    }
  }
}

SCENARIO("Test switch statement") {
  Circuit test(1);
  test.add_op<unsigned>(OpType::Ry, 1.95, {0});
  GIVEN("A circuit and a CnRy(Pi/2)") {
    Circuit circ;
    const double p = 0.5;
    WHEN("Vertex with no edges") {
      Op_ptr cnry = get_op_ptr(OpType::CnRy, p);
      circ.add_vertex(cnry);
      REQUIRE_THROWS(Transforms::decomp_controlled_Rys().apply(circ));
    }
    WHEN("Vertex with 1 edge") {
      circ.add_blank_wires(1);
      circ.add_op<unsigned>(
          OpType::CnRy, p, {0});  // automatically converted to Ry
      REQUIRE(!Transforms::decomp_controlled_Rys().apply(circ));
      REQUIRE(circ.n_vertices() == 3);  // 1 in, 1 out, 1 Ry
      REQUIRE(circ.n_gates() == 1);
      REQUIRE(circ.count_gates(OpType::Ry) == 1);
      VertexSet ry_set = circ.get_gates_of_type(OpType::Ry);
      Vertex ry = *ry_set.begin();
      REQUIRE(test_equiv_val(
          (circ.get_Op_ptr_from_Vertex(ry))->get_params().at(0), p, 4));
      REQUIRE(verify_n_qubits_for_ops(circ));
    }
    WHEN("Vertex with 2 edges") {
      circ.add_blank_wires(2);
      circ.add_op<unsigned>(OpType::CnRy, p, {0, 1});
      REQUIRE(Transforms::decomp_controlled_Rys().apply(circ));
      REQUIRE(circ.n_vertices() == 8);
      REQUIRE(circ.n_gates() == 4);
      REQUIRE(circ.count_gates(OpType::CX) == 2);
      REQUIRE(circ.count_gates(OpType::Ry) == 2);
      VertexSet ry_set = circ.get_gates_of_type(OpType::Ry);
      for (const Vertex& v : ry_set) {
        Expr param = (circ.get_Op_ptr_from_Vertex(v))->get_params().at(0);
        REQUIRE(
            (test_equiv_val(param, p / 2) || test_equiv_val(param, -p / 2)));
      }
      REQUIRE(verify_n_qubits_for_ops(circ));
    }
    WHEN("Vertex with 3 edges") {
      circ.add_blank_wires(3);
      circ.add_op<unsigned>(OpType::CnRy, p, {0, 1, 2});
      REQUIRE(Transforms::decomp_controlled_Rys().apply(circ));
      REQUIRE(circ.n_gates() == 14);
      REQUIRE(circ.count_gates(OpType::CX) == 8);
      REQUIRE(circ.count_gates(OpType::Ry) == 6);
      REQUIRE(verify_n_qubits_for_ops(circ));
    }
  }
}

SCENARIO("Test switch statement long", "[.long]") {
  Circuit test(1);
  test.add_op<unsigned>(OpType::Ry, 1.95, {0});
  const Eigen::Matrix2cd correct_block = tket_sim::get_unitary(test);
  GIVEN("A circuit and a CnRy(Pi/2)") {
    Circuit circ;
    WHEN("N-qubit CnRy gates") {
      THEN("Test with params nonzero") {
        for (unsigned N = 4; N < 10; ++N) {
          Circuit circ(N);
          std::vector<unsigned> qbs(N);
          std::iota(qbs.begin(), qbs.end(), 0);
          std::vector<Expr> params1{1.95};
          circ.add_op<unsigned>(OpType::CnRy, params1, qbs);
          REQUIRE(Transforms::decomp_controlled_Rys().apply(circ));
          const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
          const Eigen::Matrix2cd m_block =
              m.block(m.cols() - 2, m.rows() - 2, 2, 2);
          bool correctness = true;
          for (unsigned i = 0; i < 2; ++i) {
            for (unsigned j = 0; j < 2; ++j) {
              if (!approx_equal(m_block(i, j), correct_block(i, j))) {
                correctness = false;
              }
            }
          }
          REQUIRE(correctness);
          correctness = true;
          for (unsigned i = 0; i < m.rows() - 2; ++i) {
            for (unsigned j = 0; j < m.cols() - 2; ++j) {
              if (i == j) {
                if (!approx_equal(std::abs(m(i, j)), 1)) correctness = false;
              } else {
                if (!approx_equal(std::abs(m(i, j)), 0.)) correctness = false;
              }
            }
          }
          REQUIRE(correctness);
          REQUIRE(verify_n_qubits_for_ops(circ));
        }
      }
    }
  }
}

SCENARIO("Test incrementer using n borrowed qubits") {
  GIVEN("0 qbs") {
    Circuit inc = CircPool::incrementer_borrow_n_qubits(0);
    REQUIRE(inc.n_vertices() == 0);
  }
  GIVEN("A 1qb incrementer") {
    Circuit inc = CircPool::incrementer_borrow_n_qubits(1);
    REQUIRE(inc.n_gates() == 1);
    REQUIRE(inc.count_gates(OpType::X) == 1);
  }
  GIVEN("A 2qb incrementer") { REQUIRE(check_incrementer_borrow_n_qubits(2)); }
  GIVEN("A 3qb incrementer") { REQUIRE(check_incrementer_borrow_n_qubits(3)); }
  GIVEN("A 4qb incrementer") { REQUIRE(check_incrementer_borrow_n_qubits(4)); }
  GIVEN("A 5qb incrementer") { REQUIRE(check_incrementer_borrow_n_qubits(5)); }
  GIVEN("A n-qb incrementer, 5<n<10") {
    // tket_sim doesn't support computing a unitary from a 12 qubits circuit
    // hence we only test that the incrementer can be constructed as intended.
    for (unsigned n = 6; n < 10; ++n) {
      Circuit inc = CircPool::incrementer_borrow_n_qubits(n);
      REQUIRE(inc.n_qubits() == 2 * n);
      REQUIRE(inc.count_gates(OpType::CCX) == (n - 1) * 4);
      REQUIRE(Transforms::synthesise_tket().apply(inc));
    }
  }
}

SCENARIO("Test incrementer using 1 borrowed qubit") {
  GIVEN("Check the top incrementer is mapped correctly") {
    const unsigned k = 3;
    Circuit inc(2 * k);
    Circuit top_incrementer = CircPool::incrementer_borrow_n_qubits(k);
    std::vector<unsigned> top_qbs(2 * k);
    for (unsigned i = 0; i != k; ++i) {
      top_qbs[2 * i] = i + k;  // garbage qubits
      top_qbs[2 * i + 1] = i;  // qbs we are trying to increment
      inc.add_op<unsigned>(OpType::X, {i});
    }
    inc.append_qubits(top_incrementer, top_qbs);
    Transforms::decomp_CCX().apply(inc);
    const StateVector sv = tket_sim::get_statevector(inc);
    bool correct = true;
    for (unsigned i = 0; i < sv.size(); ++i) {
      if (i == 0)
        correct &= (std::abs(sv[i]) > EPS);
      else
        correct &= (std::abs(sv[i]) < ERR_EPS);
    }
    REQUIRE(correct);
  }
  GIVEN(
      "Check that the controlled bot incrementer is mapped correctly for "
      "odd qb no") {
    const unsigned j = 3;
    Circuit inc(2 * j);
    Circuit bottom_incrementer = CircPool::incrementer_borrow_n_qubits(j);
    std::vector<unsigned> bot_qbs(2 * j);
    for (unsigned i = 0; i != j; ++i) {
      bot_qbs[2 * i] = i;  // 0,2,4...n-1 //garbage qubits
      if (i != 0)
        bot_qbs[2 * i + 1] =
            i + j -
            1;  // 3,5...n //other qbs we are actually trying to increment
    }
    bot_qbs[1] = 2 * j - 1;  // incremented qubit 0 in incrementer is bottom one
    inc.add_op<unsigned>(OpType::X, {2 * j - 1});
    inc.append_qubits(bottom_incrementer, bot_qbs);
    Transforms::decomp_CCX().apply(inc);
    const StateVector sv = tket_sim::get_statevector(inc);
    bool correct = true;
    for (unsigned i = 0; i < sv.size(); ++i) {
      // |100000> -> |001000>
      if (i == 4)
        correct &= (std::abs(sv[i]) > EPS);
      else
        correct &= (std::abs(sv[i]) < ERR_EPS);
    }
    REQUIRE(correct);
  }
  GIVEN(
      "Check that the controlled bot incrementer is mapped correctly for "
      "even qb no") {
    const unsigned j = 4;
    const unsigned k = 3;
    const unsigned n = 6;
    Circuit inc(n + 1);
    for (unsigned i = k; i != n; ++i) {
      inc.add_op<unsigned>(OpType::X, {i});
    }
    Circuit bottom_incrementer = CircPool::incrementer_borrow_n_qubits(
        j - 1);  // insert incrementer over remaining qubits
    std::vector<unsigned> bot_qbs(2 * j - 2);
    for (unsigned i = 0; i != j - 1; ++i) {
      bot_qbs[2 * i] = i;  // 0,2,4...n-1 //garbage qubits
      if (i != 0)
        bot_qbs[2 * i + 1] =
            i + k -
            1;  // 3,5...n //other qbs we are actually trying to increment
    }
    bot_qbs[1] = n;  // incremented qubit 0 in incrementer is bottom one
    inc.append_qubits(bottom_incrementer, bot_qbs);
    Transforms::decomp_CCX().apply(inc);
    const StateVector sv = tket_sim::get_statevector(inc);
    bool correct = true;
    for (unsigned i = 0; i < sv.size(); ++i) {
      // |100000> -> |001000>
      if (i == 15)
        correct &= (std::abs(sv[i]) > EPS);
      else
        correct &= (std::abs(sv[i]) < ERR_EPS);
    }
    REQUIRE(correct);
  }
  GIVEN("A 0 qubit incrementer") {
    Circuit inc = CircPool::incrementer_borrow_1_qubit(0);
    REQUIRE(inc.n_qubits() == 1);
    REQUIRE(inc.n_vertices() == 2);
    REQUIRE(inc.n_gates() == 0);
  }
  GIVEN("A 1 qubit incrementer") {
    Circuit inc = CircPool::incrementer_borrow_1_qubit(1);
    REQUIRE(inc.n_qubits() == 2);
    REQUIRE(inc.n_vertices() == 5);
    REQUIRE(inc.n_gates() == 1);
  }
  GIVEN("A 2 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(2));
  }
  GIVEN("A 3 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(3));
  }
  GIVEN("A 4 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(4));
  }
  GIVEN("A 5 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(5));
  }
  GIVEN("A 6 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(6));
  }
  GIVEN("A 7 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(7));
  }
  GIVEN("A 8 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(8));
  }
  GIVEN("A 9 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(9));
  }
  GIVEN("A 10 qubit incrementer") {
    REQUIRE(check_incrementer_borrow_1_qubit(10));
  }
}

SCENARIO("Test linear depth incrementer") {
  GIVEN("0 qb") {
    Circuit circ = CircPool::incrementer_linear_depth(0);
    REQUIRE(circ.n_qubits() == 0);
  }
  GIVEN("1 qb") {
    REQUIRE(check_incrementer_linear_depth(1, 0));
    REQUIRE(check_incrementer_linear_depth(1, 1));
  }
  GIVEN("2 qbs") {
    REQUIRE(check_incrementer_linear_depth(2, 0));
    REQUIRE(check_incrementer_linear_depth(2, 1));
    REQUIRE(check_incrementer_linear_depth(2, 2));
    REQUIRE(check_incrementer_linear_depth(2, 3));
  }
  GIVEN("3 qbs") {
    REQUIRE(check_incrementer_linear_depth(3, 0));
    REQUIRE(check_incrementer_linear_depth(3, 1));
    REQUIRE(check_incrementer_linear_depth(3, 5));
    REQUIRE(check_incrementer_linear_depth(3, 7));
  }
  GIVEN("4 qbs") {
    REQUIRE(check_incrementer_linear_depth(4, 0));
    REQUIRE(check_incrementer_linear_depth(4, 1));
    REQUIRE(check_incrementer_linear_depth(4, 10));
    REQUIRE(check_incrementer_linear_depth(4, 15));
  }
  GIVEN("5 qbs") {
    REQUIRE(check_incrementer_linear_depth(5, 0));
    REQUIRE(check_incrementer_linear_depth(5, 1));
    REQUIRE(check_incrementer_linear_depth(5, 26));
    REQUIRE(check_incrementer_linear_depth(5, 31));
  }
  GIVEN("8 qbs") {
    REQUIRE(check_incrementer_linear_depth(8, 0));
    REQUIRE(check_incrementer_linear_depth(8, 1));
    REQUIRE(check_incrementer_linear_depth(8, 100));
    REQUIRE(check_incrementer_linear_depth(8, 255));
  }
}

SCENARIO("Test a CnX is decomposed correctly when bootstrapped", "[.long]") {
  GIVEN("Test CnX unitary for 3 to 9 controls") {
    for (unsigned n = 3; n < 10; ++n) {
      Circuit circ = CircPool::CnX_normal_decomp(n);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      REQUIRE(m.isApprox(get_CnX_matrix(n), ERR_EPS));
    }
  }
}

SCENARIO(
    "Test a CnX is decomposed correctly using the linear depth method",
    "[.long]") {
  GIVEN("Test CnX unitary for 0 to 9 controls") {
    for (unsigned n = 0; n < 10; ++n) {
      Eigen::MatrixXcd x = GateUnitaryMatrix::get_unitary(OpType::X, 1, {});
      Circuit circ = CircPool::CnU_linear_depth_decomp(n, x);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      REQUIRE(m.isApprox(get_CnX_matrix(n), ERR_EPS));
    }
  }
}

SCENARIO("Test a CnU is decomposed correctly using the linear depth method") {
  GIVEN("Test CnU unitary for n={0,1,2,3,5} controls") {
    for (unsigned i = 0; i < 100; i++) {
      Eigen::Matrix2cd U = random_unitary(2, i);
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_linear_depth_decomp(n, U);
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
      }
    }
  }
}

SCENARIO(
    "Test a CnU is decomposed correctly using the gray code method",
    "[.long]") {
  GIVEN("Test CnU unitary for n={0,1,2,3,5} controls") {
    for (unsigned i = 0; i < 100; i++) {
      Eigen::Matrix2cd U = random_unitary(2, i);
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_gray_code_decomp(n, U);
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
      }
    }
  }
}

static Eigen::MatrixXcd get_su2_matrix(
    const Expr& alpha, const Expr& theta, const Expr& beta) {
  Circuit c1 = Circuit(1);
  c1.add_op<unsigned>(OpType::Rz, beta, {0});
  c1.add_op<unsigned>(OpType::Ry, theta, {0});
  c1.add_op<unsigned>(OpType::Rz, alpha, {0});
  return tket_sim::get_unitary(c1);
}

SCENARIO("Test CnSU2_linear_decomp") {
  GIVEN("Test identity") {
    std::vector<std::vector<Expr>> rotations = {
        {0, 0, 4}, {0, 0, 0}, {1, 4, 3}, {1, 0, 7}, {1, 6, 1}, {1.5, -6, 4.5}};
    std::vector<unsigned> test_ns = {0, 1, 2, 3, 4, 5};
    for (auto n : test_ns) {
      for (auto angles : rotations) {
        Circuit circ =
            CircPool::CnSU2_linear_decomp(n, angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd U =
            get_su2_matrix(angles[0], angles[1], angles[2]);
        REQUIRE(U.isApprox(Eigen::Matrix2cd::Identity(), ERR_EPS));
        REQUIRE(circ.n_gates() == 0);
      }
    }
  }
  GIVEN("Test Y rotation") {
    std::vector<std::vector<Expr>> rotations = {
        {0, 0.377, 0}, {3, 4.2, -1}, {2, 1.1, 0}, {5, -0.5, -1}, {4, 1.2, 0}};
    std::vector<unsigned> test_ns = {0, 1, 2, 3, 4, 5};
    for (auto n : test_ns) {
      for (auto angles : rotations) {
        Circuit circ =
            CircPool::CnSU2_linear_decomp(n, angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd U =
            get_su2_matrix(angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
        // check the method detects pure ry
        if (n == 1) {
          REQUIRE(circ.n_gates() == 4);
          REQUIRE(circ.count_gates(OpType::CX) == 2);
          REQUIRE(circ.count_gates(OpType::Ry) == 2);
        } else if (n == 2) {
          REQUIRE(circ.n_gates() == 4);
          REQUIRE(circ.count_gates(OpType::CX) == 2);
          REQUIRE(circ.count_gates(OpType::CRy) == 2);
        }
      }
    }
  }
  GIVEN("Test rotations where W=AXBX") {
    std::vector<std::vector<Expr>> rotations = {
        {3.7, 0.377, -0.3}, {3.4, 4.2, -2.6}};
    std::vector<unsigned> test_ns = {0, 1, 2, 3, 4, 5};
    for (auto n : test_ns) {
      for (auto angles : rotations) {
        Circuit circ =
            CircPool::CnSU2_linear_decomp(n, angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd U =
            get_su2_matrix(angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
        if (n == 1) {
          REQUIRE(circ.n_gates() == 6);
          REQUIRE(circ.count_gates(OpType::CX) == 2);
          REQUIRE(circ.count_gates(OpType::Ry) == 2);
          REQUIRE(circ.count_gates(OpType::Rz) == 2);
        } else if (n == 2) {
          REQUIRE(circ.n_gates() == 6);
          REQUIRE(circ.count_gates(OpType::CX) == 2);
          REQUIRE(circ.count_gates(OpType::CRy) == 2);
          REQUIRE(circ.count_gates(OpType::CRz) == 2);
        }
      }
    }
  }
  GIVEN("Test symbolic rotations") {
    Sym a = SymTable::fresh_symbol("a");
    Expr ea(a);
    Sym b = SymTable::fresh_symbol("b");
    Expr eb(b);
    Sym c = SymTable::fresh_symbol("c");
    Expr ec(c);
    std::map<Sym, double, SymEngine::RCPBasicKeyLess> symbol_map = {
        {a, 0.3112}, {b, 1.178}, {c, -0.911}};
    std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
    for (auto n : test_ns) {
      Circuit circ = CircPool::CnSU2_linear_decomp(n, ea, eb, ec);
      const Eigen::MatrixXcd U =
          get_su2_matrix(symbol_map[a], symbol_map[b], symbol_map[c]);
      circ.symbol_substitution(symbol_map);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
    }
  }
  GIVEN("Test arbitrary rotations") {
    std::vector<std::vector<Expr>> rotations = {
        {3.3, 0.377, -0.11}, {1.3, 0, 0.13}};
    std::vector<unsigned> test_ns = {0, 1, 2, 3, 4, 5};
    for (auto n : test_ns) {
      for (auto angles : rotations) {
        Circuit circ =
            CircPool::CnSU2_linear_decomp(n, angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd U =
            get_su2_matrix(angles[0], angles[1], angles[2]);
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, U), ERR_EPS));
      }
    }
  }
}

SCENARIO(
    "Test controlled rotation gates are decomposed correctly using the gray "
    "code method",
    "[.long]") {
  GIVEN("Test CnRy for n={0,1,2,3,5} controls") {
    const Eigen::Matrix2cd ry = Gate(OpType::Ry, {Expr(3.1)}, 1).get_unitary();
    for (unsigned i = 0; i < 100; i++) {
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_gray_code_decomp(
            n, as_gate_ptr(get_op_ptr(OpType::Ry, 3.1)));
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, ry), ERR_EPS));
      }
    }
  }
  GIVEN("Test CnRx for n={0,1,2,3,5} controls") {
    const Eigen::Matrix2cd rx = Gate(OpType::Rx, {Expr(0.1)}, 1).get_unitary();
    for (unsigned i = 0; i < 100; i++) {
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_gray_code_decomp(
            n, as_gate_ptr(get_op_ptr(OpType::Rx, 0.1)));
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, rx), ERR_EPS));
      }
    }
  }
  GIVEN("Test CnRz for n={0,1,2,3,5} controls") {
    const Eigen::Matrix2cd rz = Gate(OpType::Rz, {Expr(2.7)}, 1).get_unitary();
    for (unsigned i = 0; i < 100; i++) {
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_gray_code_decomp(
            n, as_gate_ptr(get_op_ptr(OpType::Rz, 2.7)));
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, rz), ERR_EPS));
      }
    }
  }
  GIVEN("Test CnU1 for n={0,1,2,3,5} controls") {
    const Eigen::Matrix2cd u1 = Gate(OpType::U1, {Expr(1.5)}, 1).get_unitary();
    for (unsigned i = 0; i < 100; i++) {
      std::vector<unsigned> test_ns = {0, 1, 2, 3, 5};
      for (auto n : test_ns) {
        Circuit circ = CircPool::CnU_gray_code_decomp(
            n, as_gate_ptr(get_op_ptr(OpType::U1, 1.5)));
        const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
        REQUIRE(m.isApprox(get_CnU_matrix(n, u1), ERR_EPS));
      }
    }
  }
}

SCENARIO("Test a CnX is decomposed correctly using the Gray code method") {
  GIVEN("Test CnX unitary for 0 to 8 controls") {
    Circuit circ_x = CircPool::CnX_gray_decomp(0);
    REQUIRE(circ_x.n_gates() == 1);
    REQUIRE(circ_x.count_gates(OpType::X) == 1);
    Circuit circ_cx = CircPool::CnX_gray_decomp(1);
    REQUIRE(circ_cx.n_gates() == 1);
    REQUIRE(circ_cx.count_gates(OpType::CX) == 1);

    for (unsigned n = 2; n < 8; ++n) {
      Circuit circ = CircPool::CnX_gray_decomp(n);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      REQUIRE(m.isApprox(get_CnX_matrix(n), ERR_EPS));
      switch (n) {
        case 2: {
          REQUIRE(circ.count_gates(OpType::CX) <= 6);
          break;
        }
        case 3: {
          REQUIRE(circ.count_gates(OpType::CX) <= 14);
          break;
        }
        case 4: {
          REQUIRE(circ.count_gates(OpType::CX) <= 36);
          break;
        }
        case 5: {
          REQUIRE(circ.count_gates(OpType::CX) <= 92);
          break;
        }
        case 6: {
          REQUIRE(circ.count_gates(OpType::CX) <= 188);
          break;
        }
        case 7: {
          REQUIRE(circ.count_gates(OpType::CX) <= 380);
          break;
        }
      }
    }
  }
}

SCENARIO("Test decomp_arbitrary_controlled_gates") {
  GIVEN("Circuit with multi-controlled gates") {
    Circuit circ(3);

    circ.add_op<unsigned>(OpType::CnRy, 0.33, {0, 1, 2});
    circ.add_op<unsigned>(OpType::CnY, {0, 1, 2});
    circ.add_op<unsigned>(OpType::CnZ, {1, 0, 2});
    circ.add_op<unsigned>(OpType::CnX, {0, 2, 1});
    circ.add_op<unsigned>(OpType::CCX, {2, 1, 0});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::decomp_arbitrary_controlled_gates().apply(circ));
    REQUIRE(circ.count_gates(OpType::CnRy) == 0);
    REQUIRE(circ.count_gates(OpType::CnY) == 0);
    REQUIRE(circ.count_gates(OpType::CnZ) == 0);
    REQUIRE(circ.count_gates(OpType::CnX) == 0);
    REQUIRE(circ.count_gates(OpType::CCX) == 0);
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
  }

  GIVEN("Circuit without multi-controlled gates") {
    Circuit circ(3);

    circ.add_op<unsigned>(OpType::CRy, 0.33, {0, 1});
    circ.add_op<unsigned>(OpType::CRz, 0.5, {1, 2});
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    circ.add_op<unsigned>(OpType::Rx, 0.7, {0});
    REQUIRE(!Transforms::decomp_arbitrary_controlled_gates().apply(circ));
  }
}

SCENARIO("Test cnx_pairwise_decomposition") {
  GIVEN("Circuit without CnX") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CX, {1, 2});
    REQUIRE(!Transforms::cnx_pairwise_decomposition().apply(circ));
  }

  GIVEN("Circuit with C0X and C1X") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CnX, {0});
    circ.add_op<unsigned>(OpType::CnX, {0, 1});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
    REQUIRE(circ.count_gates(OpType::CnX) == 0);
  }

  GIVEN("test adding C1Z") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CnZ, {0});
    circ.add_op<unsigned>(OpType::CnZ, {0, 1});
    REQUIRE(circ.count_gates(OpType::CnZ) == 1);
  }

  GIVEN("test adding C1Y") {
    Circuit circ(2);
    circ.add_op<unsigned>(OpType::CnY, {0});
    circ.add_op<unsigned>(OpType::CnY, {0, 1});
    REQUIRE(circ.count_gates(OpType::CnY) == 1);
  }

  GIVEN("Circuit with a pair of CCX") {
    Circuit circ(3);
    circ.add_op<unsigned>(OpType::CCX, {0, 1, 2});
    circ.add_op<unsigned>(OpType::CCX, {2, 0, 1});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
    REQUIRE(Transforms::clifford_simp().apply(circ));
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
    // The CX count would normally be 12
    REQUIRE(circ.count_gates(OpType::CX) < 12);
  }

  GIVEN("CnX without overlapping qubits") {
    Circuit circ(10);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4});
    circ.add_op<unsigned>(OpType::CnX, {5, 6, 7, 8, 9});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
  }

  GIVEN("Circuit with a odd number of CnX") {
    Circuit circ(6);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4, 5});
    circ.add_op<unsigned>(OpType::CnX, {1, 2, 3, 4, 5, 0});
    circ.add_op<unsigned>(OpType::CnX, {3, 1, 4, 5, 0, 2});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
    REQUIRE(Transforms::decompose_multi_qubits_CX().apply(circ));
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
    // The CX count would normally be 240
    REQUIRE(circ.count_gates(OpType::CX) < 217);
  }

  GIVEN("Circuit with conditional CnX") {
    Circuit circ(6, 1);
    circ.add_conditional_gate<unsigned>(OpType::CnX, {}, {0, 1}, {0}, 1);
    circ.add_conditional_gate<unsigned>(
        OpType::CnX, {}, {0, 1, 2, 3, 4, 5}, {0}, 1);
    circ.add_conditional_gate<unsigned>(
        OpType::CnX, {}, {1, 2, 3, 4, 5, 0}, {0}, 1);
    circ.add_conditional_gate<unsigned>(
        OpType::CnX, {}, {3, 1, 4, 5, 0, 2}, {0}, 1);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
  }

  GIVEN("Circuit with a few more CnX") {
    Circuit circ(6);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4, 5});
    circ.add_op<unsigned>(OpType::CnX, {1, 2, 3, 4, 5, 0});
    circ.add_op<unsigned>(OpType::CnX, {2, 3, 4, 5, 0, 1});
    circ.add_op<unsigned>(OpType::CnX, {3, 4, 5, 0, 1, 2});
    circ.add_op<unsigned>(OpType::CnX, {4, 5, 0, 1, 2, 3});
    circ.add_op<unsigned>(OpType::CnX, {5, 0, 1, 2, 3, 4});
    auto u = tket_sim::get_unitary(circ);
    REQUIRE(Transforms::cnx_pairwise_decomposition().apply(circ));
    REQUIRE(Transforms::decompose_multi_qubits_CX().apply(circ));
    auto v = tket_sim::get_unitary(circ);
    REQUIRE((u - v).cwiseAbs().sum() < ERR_EPS);
    // The CX count would normally be 480
    REQUIRE(circ.count_gates(OpType::CX) < 409);
  }
}

}  // namespace test_ControlDecomp
}  // namespace tket

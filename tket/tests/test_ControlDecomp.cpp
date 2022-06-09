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

#include <catch2/catch_test_macros.hpp>
#include <numeric>

#include "Circuit/CircPool.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Transformations/ControlledGates.hpp"
#include "Transformations/OptimisationPass.hpp"
#include "Transformations/Transform.hpp"
#include "testutil.hpp"

namespace tket {
namespace test_ControlDecomp {

static bool approx_equal(const Complex& c1, const Complex& c2) {
  return (std::abs(c1 - c2) < ERR_EPS);
}

static bool check_incrementer_borrow_n_qubits(const unsigned n) {
  Circuit inc = Transforms::incrementer_borrow_n_qubits(n);
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
  Circuit inc = Transforms::incrementer_borrow_1_qubit(n);
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

SCENARIO("Test C3X and C4X decomposition") {
  GIVEN("A C3X gates") {
    Circuit circ(4);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3});
    auto u1 = tket_sim::get_unitary(circ);
    auto u2 = tket_sim::get_unitary(CircPool::C3X_normal_decomp());
    REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    Transforms::synthesise_tket().apply(circ);
    REQUIRE(circ.count_gates(OpType::CX) == 14);
  }
  GIVEN("A C4X gates") {
    Circuit circ(5);
    circ.add_op<unsigned>(OpType::CnX, {0, 1, 2, 3, 4});
    auto u1 = tket_sim::get_unitary(circ);
    auto u2 = tket_sim::get_unitary(CircPool::C4X_normal_decomp());
    REQUIRE((u1 - u2).cwiseAbs().sum() < ERR_EPS);
    Transforms::synthesise_tket().apply(circ);
    REQUIRE(circ.count_gates(OpType::CX) == 36);
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
  const Eigen::Matrix2cd correct_block = tket_sim::get_unitary(test);
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
    const double p = 0.5;
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
    Circuit inc = Transforms::incrementer_borrow_n_qubits(0);
    REQUIRE(inc.n_vertices() == 0);
  }
  GIVEN("A 1qb incrementer") {
    Circuit inc = Transforms::incrementer_borrow_n_qubits(1);
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
      Circuit inc = Transforms::incrementer_borrow_n_qubits(n);
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
    Circuit top_incrementer = Transforms::incrementer_borrow_n_qubits(k);
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
    Circuit bottom_incrementer = Transforms::incrementer_borrow_n_qubits(j);
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
    Circuit bottom_incrementer = Transforms::incrementer_borrow_n_qubits(
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
    Circuit inc = Transforms::incrementer_borrow_1_qubit(0);
    REQUIRE(inc.n_qubits() == 1);
    REQUIRE(inc.n_vertices() == 2);
    REQUIRE(inc.n_gates() == 0);
  }
  GIVEN("A 1 qubit incrementer") {
    Circuit inc = Transforms::incrementer_borrow_1_qubit(1);
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

SCENARIO("Test a CnX is decomposed correctly when bootstrapped", "[.long]") {
  GIVEN("Test CnX unitary for 3 to 9 controls") {
    for (unsigned n = 3; n < 10; ++n) {
      Circuit circ = Transforms::cnx_normal_decomp(n);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      unsigned m_size = pow(2, n + 1);
      Eigen::MatrixXcd correct_matrix =
          Eigen::MatrixXcd::Identity(m_size, m_size);
      correct_matrix(m_size - 2, m_size - 1) = 1;
      correct_matrix(m_size - 1, m_size - 2) = 1;
      correct_matrix(m_size - 2, m_size - 2) = 0;
      correct_matrix(m_size - 1, m_size - 1) = 0;
      REQUIRE(m.isApprox(correct_matrix, ERR_EPS));
    }
  }
}

SCENARIO("Test a CnX is decomposed correctly using the Gray code method") {
  GIVEN("Test CnX unitary for 0 to 8 controls") {
    Circuit circ_x = Transforms::cnx_gray_decomp(0);
    REQUIRE(circ_x.n_gates() == 1);
    REQUIRE(circ_x.count_gates(OpType::X) == 1);
    Circuit circ_cx = Transforms::cnx_gray_decomp(1);
    REQUIRE(circ_cx.n_gates() == 1);
    REQUIRE(circ_cx.count_gates(OpType::CX) == 1);

    for (unsigned n = 2; n < 8; ++n) {
      Circuit circ = Transforms::cnx_gray_decomp(n);
      const Eigen::MatrixXcd m = tket_sim::get_unitary(circ);
      unsigned m_size = pow(2, n + 1);
      Eigen::MatrixXcd correct_matrix =
          Eigen::MatrixXcd::Identity(m_size, m_size);
      correct_matrix(m_size - 2, m_size - 1) = 1;
      correct_matrix(m_size - 1, m_size - 2) = 1;
      correct_matrix(m_size - 2, m_size - 2) = 0;
      correct_matrix(m_size - 1, m_size - 1) = 0;
      REQUIRE(m.isApprox(correct_matrix, ERR_EPS));
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

}  // namespace test_ControlDecomp
}  // namespace tket

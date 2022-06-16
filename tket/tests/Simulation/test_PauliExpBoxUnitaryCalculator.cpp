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
#include <stack>

#include "Circuit/Boxes.hpp"
#include "Gate/GateUnitaryMatrix.hpp"
#include "Simulation/PauliExpBoxUnitaryCalculator.hpp"

// This is for testing PauliExpBoxUnitaryCalculator, an internal
// component of the Simulation code.

namespace tket {

namespace test_PauliExpBoxUnitaryCalculator {
struct GroupPropertyPauliTester {
  std::array<double, 3> phases;

  std::array<SparseMatrixXcd, 3> sparse_matrices;

  Eigen::MatrixXcd dense_matr1, dense_matr2;

  GroupPropertyPauliTester() {
    phases[0] = 0.123456789;
    phases[1] = 0.777666555;
    phases[2] = phases[0] + phases[1];
  }

  void test(const std::vector<Pauli>& pauli_string) {
    const auto matrix_size = get_matrix_size(pauli_string.size());

    for (unsigned ii = 0; ii < phases.size(); ++ii) {
      PauliExpBox pe_box(pauli_string, phases[ii]);
      const auto triplets = tket_sim::internal::get_triplets(pe_box);
      sparse_matrices[ii] = get_sparse_square_matrix(triplets, matrix_size);
    }

    // Check group property (homomorphism).
    dense_matr1 = sparse_matrices[0] * sparse_matrices[1];
    dense_matr2 = sparse_matrices[2];
    CHECK(dense_matr1.isApprox(dense_matr2));

    // Commutation.
    dense_matr1 = sparse_matrices[1] * sparse_matrices[0];
    CHECK(dense_matr1.isApprox(dense_matr2));

    // Inverse.
    PauliExpBox pe_box(pauli_string, -phases[0]);
    const auto triplets = tket_sim::internal::get_triplets(pe_box);
    sparse_matrices[1] = get_sparse_square_matrix(triplets, matrix_size);
    dense_matr1 = sparse_matrices[0] * sparse_matrices[1];
    CHECK(dense_matr1.isApprox(
        Eigen::MatrixXcd::Identity(matrix_size, matrix_size)));
  }
};

// Directly calculate the tensor product with Eigen
// to recalculate the value.
struct DirectTensorProductTester {
  const double phase;
  const double cc;
  const double ss;
  std::map<Pauli, Eigen::MatrixXcd> pauli_matrices;
  SparseMatrixXcd sparse_matrix;
  Eigen::MatrixXcd dense_matrix;
  Eigen::MatrixXcd recalculated_result;

  DirectTensorProductTester()
      : phase(-1.234567),
        cc(std::cos(0.5 * PI * phase)),
        ss(std::sin(0.5 * PI * phase)) {
    auto const_2x2_matrix = [](Complex tl, Complex tr, Complex bl, Complex br) {
      CmplxSpMat m(2, 2);
      if (tl != czero) {
        m.insert(0, 0) = tl;
      }
      if (tr != czero) {
        m.insert(0, 1) = tr;
      }
      if (bl != czero) {
        m.insert(1, 0) = bl;
      }
      if (br != czero) {
        m.insert(1, 1) = br;
      }
      return m;
    };
    pauli_matrices[Pauli::I] = const_2x2_matrix(1, 0, 0, 1);
    pauli_matrices[Pauli::X] = const_2x2_matrix(0, 1, 1, 0);
    pauli_matrices[Pauli::Y] = const_2x2_matrix(0, -i_, i_, 0);
    pauli_matrices[Pauli::Z] = const_2x2_matrix(1, 0, 0, -1);
  }

  void test(const std::vector<Pauli>& pauli_string) {
    const auto matrix_size = get_matrix_size(pauli_string.size());
    PauliExpBox pe_box(pauli_string, phase);
    const auto triplets = tket_sim::internal::get_triplets(pe_box);
    sparse_matrix = get_sparse_square_matrix(triplets, matrix_size);
    dense_matrix = sparse_matrix;

    recalculated_result.resize(1, 1);
    recalculated_result(0, 0) = 1;
    for (auto pauli : pauli_string) {
      // Apparently, we need eval() to workaround an aliasing bug.
      recalculated_result =
          Eigen::kroneckerProduct(recalculated_result, pauli_matrices.at(pauli))
              .eval();
    }
    const Eigen::MatrixXcd identity =
        Eigen::MatrixXcd::Identity(matrix_size, matrix_size);
    CHECK(identity.isApprox(recalculated_result * recalculated_result));
    recalculated_result = cc * identity - (i_ * ss) * recalculated_result;
    CHECK(recalculated_result.isApprox(dense_matrix));
  }
};

SCENARIO("Check all Pauli strings of length <= 4") {
  GroupPropertyPauliTester group_property_tester;
  DirectTensorProductTester direct_product_tester;

  std::stack<std::vector<Pauli>> pauli_strings;
  pauli_strings.push({});
  const std::array<Pauli, 3> non_identity_paulis{Pauli::X, Pauli::Y, Pauli::Z};

  unsigned count = 0;
  while (!pauli_strings.empty()) {
    group_property_tester.test(pauli_strings.top());
    direct_product_tester.test(pauli_strings.top());
    ++count;
    if (pauli_strings.top().size() >= 4) {
      pauli_strings.pop();
    } else {
      pauli_strings.top().push_back(Pauli::I);
      auto copy = pauli_strings.top();
      for (auto pauli : non_identity_paulis) {
        copy.back() = pauli;
        pauli_strings.push(copy);
      }
    }
  }
  // count the tests! Did we REALLY test every expected string?
  CHECK(count == 1 + 4 + 16 + 64 + 256);
}

// Use other OpTypes which give equivalent unitaries to some Pauli Exp Boxes.
struct CompareWithSimulatorPauliTester {
  const std::vector<double> phase;
  SparseMatrixXcd result;
  Eigen::MatrixXcd dense_result;
  Eigen::MatrixXcd simulator_result;

  CompareWithSimulatorPauliTester() : phase({0.987654321}) {}

  // Find an OpType which gives a gate with equivalent unitary to
  // the Pauli Exp Box.
  void test(OpType type, const std::vector<Pauli>& paulis) {
    const auto matr_size = get_matrix_size(paulis.size());

    // For convenience, treat "noop" as a special marker:
    // there's no single OpType to give a multiple of the identity.
    if (type == OpType::noop) {
      simulator_result = std::polar(1.0, -0.5 * PI * phase[0]) *
                         Eigen::MatrixXcd::Identity(matr_size, matr_size);
    } else {
      simulator_result =
          GateUnitaryMatrix::get_unitary(type, paulis.size(), phase);
    }
    REQUIRE(simulator_result.rows() == matr_size);
    REQUIRE(simulator_result.cols() == matr_size);

    PauliExpBox pe_box(paulis, phase[0]);
    const auto triplets = tket_sim::internal::get_triplets(pe_box);
    result = get_sparse_square_matrix(triplets, matr_size);
    dense_result = result;
    CHECK(dense_result.isApprox(simulator_result));
  }
};

SCENARIO("Check some length <= 2 Pauli strings using equivalent gates") {
  CompareWithSimulatorPauliTester tester;

  // Pauli -> equivalent 1 and 2 qubit gates.
  const std::map<Pauli, std::pair<OpType, OpType>> data{
      // Using noop as a special case
      {Pauli::I, {OpType::noop, OpType::noop}},
      {Pauli::X, {OpType::Rx, OpType::XXPhase}},
      {Pauli::Y, {OpType::Ry, OpType::YYPhase}},
      {Pauli::Z, {OpType::Rz, OpType::ZZPhase}}};
  std::vector<Pauli> paulis;

  for (const auto& entry : data) {
    // II, XX, YY, ZZ
    paulis.resize(2);
    paulis[0] = paulis[1] = entry.first;
    tester.test(entry.second.second, paulis);

    // Single Paulis
    paulis.resize(1);
    tester.test(entry.second.first, paulis);
  }
}

}  // namespace test_PauliExpBoxUnitaryCalculator
}  // namespace tket

// Copyright 2019-2024 Cambridge Quantum Computing
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

#include "tket/Circuit/Simulation/CircuitSimulator.hpp"
#include "tket/Diagonalisation/Diagonalisation.hpp"

namespace tket {

namespace test_Diagonalisation {

SCENARIO("Matrix tests for reducing a Pauli to Z") {
  std::list<Pauli> test_paulis{Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
  std::list<CXConfigType> test_configs = {
      CXConfigType::Snake, CXConfigType::Tree, CXConfigType::Star,
      CXConfigType::MultiQGate};
  for (const Pauli& p : test_paulis) {
    // If p is Pauli::I, it is dropped from the sparse representation in the
    // constructor, so need to add Qubit(3) to the circuit and sparse matrix
    // explicitly
    SpPauliStabiliser pt({Pauli::X, Pauli::Y, Pauli::Z, p});
    CmplxSpMat pt_u = pt.to_sparse_matrix(4);
    Eigen::MatrixXcd pt_ud = pt_u;
    for (const CXConfigType& config : test_configs) {
      std::pair<Circuit, Qubit> diag = reduce_pauli_to_z(pt, config);
      if (p == Pauli::I) diag.first.add_qubit(Qubit(3));
      Eigen::MatrixXcd diag_u = tket_sim::get_unitary(diag.first);
      CmplxSpMat z_u = SpPauliString(diag.second, Pauli::Z).to_sparse_matrix(4);
      Eigen::MatrixXcd z_ud = z_u;
      CHECK((z_ud * diag_u * pt_ud).isApprox(diag_u));
    }
  }
}

SCENARIO("Matrix tests for reducing two anticommuting Paulis to Z X") {
  std::list<Pauli> non_trivials{Pauli::X, Pauli::Y, Pauli::Z};
  std::list<CXConfigType> test_configs = {
      CXConfigType::Snake, CXConfigType::Tree, CXConfigType::Star,
      CXConfigType::MultiQGate};
  // Loop through all commuting options for two qubits
  for (const Pauli& p0 : non_trivials) {
    for (const Pauli& p1 : non_trivials) {
      SpPauliStabiliser p({Pauli::Z, p0, p1, Pauli::Z});
      CmplxSpMat p_u = p.to_sparse_matrix();
      Eigen::MatrixXcd p_ud = p_u;
      for (const Pauli& q0 : non_trivials) {
        for (const Pauli& q1 : non_trivials) {
          SpPauliStabiliser q({Pauli::X, q0, q1, Pauli::Z});
          if (p.commutes_with(q)) continue;
          CmplxSpMat q_u = q.to_sparse_matrix();
          Eigen::MatrixXcd q_ud = q_u;
          for (const CXConfigType& config : test_configs) {
            std::pair<Circuit, Qubit> diag =
                reduce_anticommuting_paulis_to_z_x(p, q, config);
            Eigen::MatrixXcd diag_u = tket_sim::get_unitary(diag.first);
            CmplxSpMat z_u =
                SpPauliString(diag.second, Pauli::Z).to_sparse_matrix(4);
            Eigen::MatrixXcd z_ud = z_u;
            CmplxSpMat x_u =
                SpPauliString(diag.second, Pauli::X).to_sparse_matrix(4);
            Eigen::MatrixXcd x_ud = x_u;
            CHECK((z_ud * diag_u * p_ud).isApprox(diag_u));
            CHECK((x_ud * diag_u * q_ud).isApprox(diag_u));
          }
        }
      }
    }
  }
}

SCENARIO("Matrix tests for reducing two commuting Paulis to Z X") {
  std::list<Pauli> paulis{Pauli::I, Pauli::X, Pauli::Y, Pauli::Z};
  std::list<CXConfigType> test_configs = {
      CXConfigType::Snake, CXConfigType::Tree, CXConfigType::Star,
      CXConfigType::MultiQGate};
  // Loop through all commuting options for two qubits
  for (const Pauli& p0 : paulis) {
    for (const Pauli& p1 : paulis) {
      SpPauliStabiliser p({Pauli::Z, p0, p1, Pauli::Z});
      CmplxSpMat p_u = p.to_sparse_matrix(4);
      Eigen::MatrixXcd p_ud = p_u;
      for (const Pauli& q0 : paulis) {
        for (const Pauli& q1 : paulis) {
          SpPauliStabiliser q({Pauli::Z, q0, q1, Pauli::I});
          if (!p.commutes_with(q)) continue;
          CmplxSpMat q_u = q.to_sparse_matrix(4);
          Eigen::MatrixXcd q_ud = q_u;
          for (const CXConfigType& config : test_configs) {
            std::tuple<Circuit, Qubit, Qubit> diag =
                reduce_commuting_paulis_to_zi_iz(p, q, config);
            Circuit& circ = std::get<0>(diag);
            // In cases with matching Pauli::Is, the circuit produced may not
            // include all qubits
            for (unsigned i = 0; i < 4; ++i) circ.add_qubit(Qubit(i), false);
            Eigen::MatrixXcd diag_u = tket_sim::get_unitary(circ);
            CmplxSpMat zi_u =
                SpPauliString(std::get<1>(diag), Pauli::Z).to_sparse_matrix(4);
            Eigen::MatrixXcd zi_ud = zi_u;
            CmplxSpMat iz_u =
                SpPauliString(std::get<2>(diag), Pauli::Z).to_sparse_matrix(4);
            Eigen::MatrixXcd iz_ud = iz_u;
            CHECK((zi_ud * diag_u * p_ud).isApprox(diag_u));
            CHECK((iz_ud * diag_u * q_ud).isApprox(diag_u));
          }
        }
      }
    }
  }
}

}  // namespace test_Diagonalisation
}  // namespace tket
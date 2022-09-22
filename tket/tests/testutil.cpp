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

#include "testutil.hpp"

#include <catch2/catch_test_macros.hpp>

#include "Circuit/Circuit.hpp"
#include "Simulation/CircuitSimulator.hpp"
#include "Simulation/ComparisonFunctions.hpp"
#include "Utils/EigenConfig.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {

bool test_statevector_comparison(
    const Circuit& circ1, const Circuit& circ2, bool projective) {
  const StateVector s1 = tket_sim::get_statevector(circ1);
  const StateVector s2 = tket_sim::get_statevector(circ2);
  return tket_sim::compare_statevectors_or_unitaries(
      s1, s2,
      projective ? tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE
                 : tket_sim::MatrixEquivalence::EQUAL);
}

bool test_unitary_comparison(
    const Circuit& circ1, const Circuit& circ2, bool projective) {
  const Eigen::MatrixXcd m1 = tket_sim::get_unitary(circ1);
  const Eigen::MatrixXcd m2 = tket_sim::get_unitary(circ2);
  return tket_sim::compare_statevectors_or_unitaries(
      m1, m2,
      projective ? tket_sim::MatrixEquivalence::EQUAL_UP_TO_GLOBAL_PHASE
                 : tket_sim::MatrixEquivalence::EQUAL);
}

bool verify_n_qubits_for_ops(const Circuit& circ) {
  for (const Command& com : circ) {
    if (com.get_op_ptr()->n_qubits() != com.get_args().size()) {
      return false;
    }
  }
  return true;
}

void add_2qb_gates(
    Circuit& circ, OpType op_type,
    const std::vector<std::pair<unsigned, unsigned>>& qubit_pairs) {
  std::vector<unsigned> qbs(2);
  for (const auto& pair : qubit_pairs) {
    qbs[0] = pair.first;
    qbs[1] = pair.second;
    circ.add_op<unsigned>(op_type, qbs);
  }
}

void add_1qb_gates(
    Circuit& circ, OpType op_type, const std::vector<unsigned>& qubits) {
  std::vector<unsigned> qbs(1);
  for (unsigned qubit : qubits) {
    qbs[0] = qubit;
    circ.add_op<unsigned>(op_type, qbs);
  }
}

void reassign_boundary(Circuit& circ_, std::optional<node_vector_t> nodes) {
  qubit_vector_t allqs = circ_.all_qubits();
  unit_map_t rename_map;
  for (const Qubit& q : allqs) {
    Node n;
    if (nodes.has_value()) {
      n = nodes->at(q.index().front());
    } else {
      n = Node(q.index().front());
    }
    rename_map.insert({q, n});
  }
  circ_.rename_units(rename_map);
}

void check_command_types(
    const Circuit& circ, const std::vector<OpType>& expected_types) {
  const std::vector<Command> coms = circ.get_commands();
  REQUIRE(coms.size() == expected_types.size());
  for (unsigned nn = 0; nn < coms.size(); ++nn) {
    INFO("circuit " << circ << ", nn=" << nn);
    REQUIRE_NOTHROW(coms[nn].to_str());
    REQUIRE(coms[nn].get_op_ptr()->get_type() == expected_types[nn]);
  }
}

Eigen::MatrixXcd random_unitary(unsigned n, int seed) {
  std::srand(seed);
  Eigen::MatrixXcd A = Eigen::MatrixXcd::Random(n, n);
  Eigen::MatrixXcd H = A + A.adjoint();
  // H is Hermitian, so exp(iH) is unitary.
  return (i_ * H).exp();
}

}  // namespace tket

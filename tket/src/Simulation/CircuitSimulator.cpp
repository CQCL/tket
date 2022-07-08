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

#include "CircuitSimulator.hpp"

#include <sstream>

#include "Circuit/Circuit.hpp"
#include "DecomposeCircuit.hpp"
#include "Gate/GateUnitaryMatrixError.hpp"
#include "GateNodesBuffer.hpp"
#include "Utils/Expression.hpp"
#include "Utils/UnitID.hpp"

namespace tket {
namespace tket_sim {

Eigen::MatrixXcd get_unitary(
    const Circuit& circ, double abs_epsilon, unsigned max_number_of_qubits) {
  const auto matr_size = get_matrix_size(circ.n_qubits());
  Eigen::MatrixXcd result = Eigen::MatrixXcd::Identity(matr_size, matr_size);
  apply_unitary(circ, result, abs_epsilon, max_number_of_qubits);
  return result;
}

static void apply_unitary_may_throw(
    const Circuit& circ, Eigen::MatrixXcd& matr, double abs_epsilon,
    unsigned max_number_of_qubits) {
  if (circ.n_qubits() > max_number_of_qubits) {
    throw GateUnitaryMatrixError(
        "Circuit to simulate has too many qubits",
        GateUnitaryMatrixError::Cause::TOO_MANY_QUBITS);
  }
  if (matr.cols() <= 0) {
    throw GateUnitaryMatrixError(
        "M has no columns", GateUnitaryMatrixError::Cause::INPUT_ERROR);
  }
  const auto full_matr_size = get_matrix_size(circ.n_qubits());
  if (matr.rows() != full_matr_size) {
    throw GateUnitaryMatrixError(
        "M has wrong number of rows",
        GateUnitaryMatrixError::Cause::INPUT_ERROR);
  }
  internal::GateNodesBuffer buffer(matr, abs_epsilon);
  internal::decompose_circuit(circ, buffer, abs_epsilon);

  // Apply qubit permutation:
  const qubit_map_t perm = circ.implicit_qubit_permutation();
  std::map<Qubit, unsigned> q_indices;
  std::map<unsigned, unsigned> uq_map;
  unsigned qi = 0;
  for (const std::pair<const Qubit, Qubit>& pair : perm) {
    q_indices.insert({pair.first, qi});
    ++qi;
  }
  for (const std::pair<const Qubit, Qubit>& pair : perm) {
    unsigned in = q_indices[pair.first];
    unsigned out = q_indices[pair.second];
    uq_map.insert({in, out});
  }
  matr = lift_perm(uq_map) * matr;
}

void apply_unitary(
    const Circuit& circ, Eigen::MatrixXcd& matr, double abs_epsilon,
    unsigned max_number_of_qubits) {
  try {
    apply_unitary_may_throw(circ, matr, abs_epsilon, max_number_of_qubits);
  } catch (const GateUnitaryMatrixError& e) {
    const auto full_matr_size = get_matrix_size(circ.n_qubits());
    std::stringstream ss;
    ss << "Error trying to simulate circuit " << circ << " with "
       << circ.n_qubits() << " qubits, " << circ.get_commands().size()
       << " commands; U is size " << full_matr_size << "x" << full_matr_size
       << ", premultiplying M with " << matr.rows() << " rows, " << matr.cols()
       << " cols: " << e.what();
    if (e.cause == GateUnitaryMatrixError::Cause::GATE_NOT_IMPLEMENTED) {
      throw CircuitInvalidity(ss.str());
    } else if (e.cause == GateUnitaryMatrixError::Cause::SYMBOLIC_PARAMETERS) {
      throw SymbolsNotSupported(ss.str());
    } else {
      throw e;
    }
  }
}

Eigen::VectorXcd get_statevector(
    const Circuit& circ, double abs_epsilon, unsigned max_number_of_qubits) {
  Eigen::MatrixXcd result =
      Eigen::MatrixXcd::Zero(get_matrix_size(circ.n_qubits()), 1);
  result(0, 0) = 1.0;
  apply_unitary(circ, result, abs_epsilon, max_number_of_qubits);
  return result;
}

}  // namespace tket_sim
}  // namespace tket

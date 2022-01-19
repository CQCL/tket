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

#include "GateUnitaryMatrixUtils.hpp"

#include <sstream>

#include "Gate/Gate.hpp"
#include "GateUnitaryMatrixError.hpp"

namespace tket {
namespace internal {

Eigen::Matrix4cd GateUnitaryMatrixUtils::get_controlled_gate_unitary(
    const Eigen::Matrix2cd& u) {
  Eigen::Matrix4cd matr = Eigen::Matrix4cd::Identity();
  matr.block(2, 2, 2, 2) = u;
  return matr;
}

Eigen::MatrixXcd
GateUnitaryMatrixUtils::get_multi_controlled_gate_dense_unitary(
    const Eigen::MatrixXcd& u, unsigned int number_of_qubits) {
  unsigned int size;
  try {
    size = get_matrix_size(number_of_qubits);
  } catch (const std::exception& e) {
    throw GateUnitaryMatrixError(
        e.what(), GateUnitaryMatrixError::Cause::TOO_MANY_QUBITS);
  }
  const auto throw_with_message = [&](const std::string& message) {
    std::stringstream ss;
    ss << "multi_controlled_gate with " << number_of_qubits
       << " qubits "
          "(final matrix size "
       << size << "x" << size
       << "), "
          "for unitary matrix U with "
       << u.cols() << " cols, " << u.rows() << ": " << message;
    throw GateUnitaryMatrixError(
        ss.str(), GateUnitaryMatrixError::Cause::INPUT_ERROR);
  };
  if (u.cols() != u.rows()) {
    throw_with_message("matrix U not square");
  }
  if (u.cols() == 0) {
    throw_with_message("zero size matrix U");
  }
  if (!(number_of_qubits >= 1 && size >= 2)) {
    throw_with_message("must have at least 1 qubit");
  }
  if (!(size >= u.cols())) {
    throw_with_message("input U is too large for the final number of qubits");
  }
  // A trick: we should check that U is of size 2^k * 2^k.
  // But actually, we know the full size is 2^N * 2^N for N >= k,
  // so we just check for factors.
  const unsigned int ratio = size / u.cols();
  if (ratio * u.cols() != size) {
    std::stringstream ss;
    ss << "input U number of columns is not a power of 2 (" << u.cols()
       << " doesn't divide " << size << ")";
    throw_with_message(ss.str());
  }
  Eigen::MatrixXcd matr = Eigen::MatrixXcd::Identity(size, size);
  matr.block(size - u.cols(), size - u.cols(), u.cols(), u.cols()) = u;
  return matr;
}

std::string GateUnitaryMatrixUtils::get_error_prefix(
    const std::string& name, unsigned number_of_qubits,
    const std::vector<double>& parameters) {
  std::stringstream ss;
  ss << "GateUnitaryMatrix for op " << name << " acting on " << number_of_qubits
     << " qubits, taking " << parameters.size() << " parameters:\n";

  for (unsigned nn = 0; nn < parameters.size(); ++nn) {
    if (nn >= 10) {
      ss << "...";
      break;
    }
    ss << "param[" << nn << "] = " << parameters[nn] << "\n";
  }
  return ss.str();
}

std::string GateUnitaryMatrixUtils::get_error_prefix(
    OpType op_type, unsigned number_of_qubits,
    const std::vector<double>& parameters) {
  const OpDesc desc(op_type);
  return get_error_prefix(desc.name(), number_of_qubits, parameters);
}

void GateUnitaryMatrixUtils::check_and_throw_upon_wrong_number_of_parameters(
    OpType op_type, unsigned number_of_qubits,
    const std::vector<double>& parameters,
    unsigned expected_number_of_parameters) {
  if (parameters.size() == expected_number_of_parameters) {
    return;
  }
  std::stringstream ss;
  ss << get_error_prefix(op_type, number_of_qubits, parameters)
     << "wrong number of parameters (expected " << expected_number_of_parameters
     << ")";
  throw GateUnitaryMatrixError(
      ss.str(), GateUnitaryMatrixError::Cause::INPUT_ERROR);
}

std::vector<double> GateUnitaryMatrixUtils::get_checked_parameters(
    const Gate& gate) {
  const std::vector<Expr> parameter_expressions = gate.get_params();
  const unsigned int number_of_qubits = gate.n_qubits();
  std::vector<double> parameters(parameter_expressions.size());
  for (unsigned nn = 0; nn < parameters.size(); ++nn) {
    const auto optional_value = eval_expr(parameter_expressions[nn]);
    if (!optional_value) {
      std::stringstream ss;
      ss << get_error_prefix(gate.get_name(), number_of_qubits, parameters)
         << "parameter[" << nn << "] is symbolic";
      throw GateUnitaryMatrixError(
          ss.str(), GateUnitaryMatrixError::Cause::SYMBOLIC_PARAMETERS);
    }
    const double value = optional_value.value();
    if (!std::isfinite(value)) {
      std::stringstream ss;
      ss << get_error_prefix(gate.get_name(), number_of_qubits, parameters)
         << "parameter[" << nn << "] has non-finite value " << value;
      throw GateUnitaryMatrixError(
          ss.str(), GateUnitaryMatrixError::Cause::NON_FINITE_PARAMETER);
    }
    parameters[nn] = value;
  }
  return parameters;
}

}  // namespace internal
}  // namespace tket

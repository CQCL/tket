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

#pragma once

#include <string>

#include "OpType/OpType.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {
class Gate;
namespace internal {

struct GateUnitaryMatrixUtils {
  /** For an arbitrary 1-qubit gate represented by a unitary,
   *  compute the unitary of the 2-qubit gate, in ILO-BE convention,
   *  obtained by adding a control qubit before the target qubit.
   *  @param u The 2x2 unitary (although unitarity is not checked).
   */
  static Eigen::Matrix4cd get_controlled_gate_unitary(
      const Eigen::Matrix2cd& u);

  /** For an arbitrary unitary matrix representing an m-qubit gate,
   *  (although unitarity is not checked), in ILO-BE convention,
   *  return the matrix of the N-qubit gate obtained by adding control qubits
   *  before the target qubits.
   *  @param u General unitary matrix (although unitarity is not checked).
   *  @param number_of_qubits Total number of qubits (controls and targets).
   */
  static Eigen::MatrixXcd get_multi_controlled_gate_dense_unitary(
      const Eigen::MatrixXcd& u, unsigned int number_of_qubits);

  /** Used to get consistent error messages. */
  static std::string get_error_prefix(
      const std::string& op_name, unsigned number_of_qubits,
      const std::vector<double>& parameters);

  static std::string get_error_prefix(
      OpType op_type, unsigned number_of_qubits,
      const std::vector<double>& parameters);

  /** Throws if the number of parameters doesn't match expectations. */
  static void check_and_throw_upon_wrong_number_of_parameters(
      OpType op_type, unsigned number_of_qubits,
      const std::vector<double>& parameters,
      unsigned expected_number_of_parameters);

  /** Converts the parameters into doubles, with checks that they all
   *  have finite numerical values.
   *  Throws a GateUnitaryMatrixError upon error.
   */
  static std::vector<double> get_checked_parameters(const Gate& gate);
};

}  // namespace internal
}  // namespace tket

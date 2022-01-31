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

#include "OpType/OpType.hpp"
#include "Utils/MatrixAnalysis.hpp"

namespace tket {
namespace internal {

/** Gates taking a variable number of qubits are treated as a special case.
 *  This class checks the type, and knows how to get the unitary.
 */
class GateUnitaryMatrixVariableQubits {
 public:
  explicit GateUnitaryMatrixVariableQubits(OpType op_type);

  bool is_known_type() const;
  unsigned get_number_of_parameters() const;

  /** Call this only if known_type() returned true,
   *  and if the number of parameters matches the parameters size.
   *  Returns the unitary for the type passed into the constructor.
   *  Uses ILO-BE convention.
   *  Does not check for finite values.
   */
  Eigen::MatrixXcd get_dense_unitary(
      unsigned number_of_qubits, const std::vector<double>& parameters) const;

 private:
  const OpType op_type;
  bool known_type;
  unsigned number_of_parameters;
};

}  // namespace internal
}  // namespace tket

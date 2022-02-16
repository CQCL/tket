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

#include "GateUnitaryMatrixVariableQubits.hpp"

#include "GateUnitaryMatrixImplementations.hpp"
#include "Utils/Assert.hpp"

namespace tket {
namespace internal {

GateUnitaryMatrixVariableQubits::GateUnitaryMatrixVariableQubits(
    OpType op_type_)
    : op_type(op_type_), known_type(true), number_of_parameters(0) {
  switch (op_type) {
    case OpType::CnRy:
      // Fall through.
    case OpType::PhaseGadget:
      number_of_parameters = 1;
      break;
    case OpType::CnX:
      break;
    case OpType::NPhasedX:
      number_of_parameters = 2;
      break;
    default:
      known_type = false;
  }
}

bool GateUnitaryMatrixVariableQubits::is_known_type() const {
  return known_type;
}

unsigned GateUnitaryMatrixVariableQubits::get_number_of_parameters() const {
  return number_of_parameters;
}

Eigen::MatrixXcd GateUnitaryMatrixVariableQubits::get_dense_unitary(
    unsigned number_of_qubits, const std::vector<double>& parameters) const {
  // This class is internal only, so an assert is OK.
  TKET_ASSERT(known_type);
  TKET_ASSERT(parameters.size() == number_of_parameters);
  switch (parameters.size()) {
    case 0:
      TKET_ASSERT(op_type == OpType::CnX);
      return GateUnitaryMatrixImplementations::CnX(number_of_qubits);
    case 1:
      if (op_type == OpType::CnRy) {
        return GateUnitaryMatrixImplementations::CnRy(
            number_of_qubits, parameters[0]);
      } else {
        TKET_ASSERT(op_type == OpType::PhaseGadget);
        return GateUnitaryMatrixImplementations::PhaseGadget(
            number_of_qubits, parameters[0]);
      }
    case 2:
      TKET_ASSERT(op_type == OpType::NPhasedX);
      return GateUnitaryMatrixImplementations::NPhasedX(
          number_of_qubits, parameters[0], parameters[1]);
    default:
      TKET_ASSERT(false);
  }
}

}  // namespace internal
}  // namespace tket

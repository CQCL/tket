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

#include <stdexcept>
#include <string>

namespace tket {

/** We do not wish to couple tightly to exception classes elsewhere.
 *  Let the caller catch this and decide what to do,
 *  based upon the cause of the error.
 */
struct GateUnitaryMatrixError : public std::runtime_error {
  enum class Cause {
    SYMBOLIC_PARAMETERS,
    GATE_NOT_IMPLEMENTED,
    // The matrix would be too big
    TOO_MANY_QUBITS,
    NON_FINITE_PARAMETER,
    // The routines were passed invalid input
    // (wrong number of qubits/arguments/etc.) somehow.
    INPUT_ERROR
  };

  const Cause cause;

  explicit GateUnitaryMatrixError(const std::string& message, Cause cause);
};

}  // namespace tket

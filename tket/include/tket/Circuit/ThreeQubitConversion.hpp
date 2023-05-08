// Copyright 2019-2023 Cambridge Quantum Computing
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

#include "Circuit.hpp"
#include "Utils/EigenConfig.hpp"

namespace tket {

/**
 * Synthesise a 3-qubit circuit from an arbitrary 8x8 unitary.
 *
 * The returned circuit consists of CX, TK1, H, Ry and Rz gates only. It
 * contains a maximum of 20 CX gates.
 *
 * @param U unitary matrix in \ref BasisOrder::ilo
 *
 * @return circuit implementing the unitary
 */
Circuit three_qubit_synthesis(const Eigen::MatrixXcd &U);

/**
 * Synthesise a 3-qubit circuit from an arbitrary 8x8 unitary.
 *
 * The returned circuit consists of TK2 and 1-qubit gates only. It contains a
 * maximum of 15 TK2 gates.
 *
 * @param U unitary matrix in \ref BasisOrder::ilo
 *
 * @return circuit implementing the unitary
 */
Circuit three_qubit_tk_synthesis(const Eigen::MatrixXcd &U);

/**
 * Convert a 3-qubit circuit to its corresponding unitary matrix.
 *
 * @param c pure quantum circuit with 3 qubits
 *
 * @pre \p c is composed of 1- and 2-qubit gates only
 * @pre \p c has no symbolic parameters
 *
 * @return 8x8 unitary matrix in \ref BasisOrder::ilo
 */
Eigen::MatrixXcd get_3q_unitary(const Circuit &c);

}  // namespace tket

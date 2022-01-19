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

#include "Circuit/Circuit.hpp"
#include "Utils/PauliStrings.hpp"

namespace tket {

/**
 * Append a Pauli gadget to the end of a given circuit.
 * Automatically uses Snake CX configuration
 *
 * @param circ circuit to append to
 * @param pauli Pauli operators and their respective qubits
 * @param angle angle in half-turns
 * @param cx_config which type of CX configuration to decompose into
 */
void append_single_pauli_gadget(
    Circuit& circ, const QubitPauliTensor& pauli, Expr angle,
    CXConfigType cx_config = CXConfigType::Snake);

/**
 * Append a pair of Pauli gadgets to the end of a given circuit.
 * (shallow) Uses an adapted arrangement of CX that gives balanced trees
 * over the matching qubits to improve depth. Better performance
 * is not guaranteed as CXs may not align for cancellation and
 * it can be harder to route.
 * (!shallow) Uses the original method with naive arrangement of CXs.
 *
 * @param circ circuit to append to
 * @param pauli0 first Pauli string
 * @param angle0 angle for \p pauli0 (half-turns)
 * @param pauli1 second Pauli string
 * @param angle1 angle for \p pauli1 (half-turns)
 * @param cx_config which type of CX configuration to decompose into
 */
void append_pauli_gadget_pair(
    Circuit& circ, QubitPauliTensor pauli0, Expr angle0,
    QubitPauliTensor pauli1, Expr angle1,
    CXConfigType cx_config = CXConfigType::Snake);

}  // namespace tket

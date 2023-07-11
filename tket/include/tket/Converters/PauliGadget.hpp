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

#include <stdexcept>

#include "tket/Architecture/Architecture.hpp"
#include "tket/Circuit/Circuit.hpp"
#include "tket/Utils/PauliStrings.hpp"

namespace tket {

class ImplicitPermutationNotAllowed : public std::logic_error {
 public:
  explicit ImplicitPermutationNotAllowed(const std::string& message)
      : std::logic_error(message) {}
};

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
 * Append a Pauli gadget to the end of a given circuit as a
 * PauliExpBox.
 * Automatically uses Snake CX configuration
 *
 * @param circ circuit to append to
 * @param pauli Pauli operators and their respective qubits
 * @param angle angle in half-turns
 * @param cx_config which type of CX configuration to decompose into
 */
void append_single_pauli_gadget_as_pauli_exp_box(
    Circuit& circ, const QubitPauliTensor& pauli, Expr angle,
    CXConfigType cx_config);

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

void append_pauli_gadget_pair_as_box(
    Circuit& circ, const QubitPauliTensor& pauli0, Expr angle0,
    const QubitPauliTensor& pauli1, Expr angle1, CXConfigType cx_config);

void append_commuting_pauli_gadget_set_as_box(
    Circuit& circ, const std::list<std::pair<QubitPauliTensor, Expr>>& gadgets,
    CXConfigType cx_config);

/**
 * @brief AAS version of append_commuting_pauli_gadget_set_as_box
 * however apprend a decomposed circuit
 *
 * @param circ
 * @param gadgets
 * @param arch
 */
void append_aased_commuting_pauli_gadget_set(
    Circuit& circ, std::list<std::pair<QubitPauliTensor, Expr>> gadgets,
    const Architecture& arch,
    const std::function<Circuit(const Architecture&, const Circuit&)>& aas_func
    );

}  // namespace tket

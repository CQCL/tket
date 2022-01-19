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

#include "Circuit/Boxes.hpp"
#include "Circuit/Circuit.hpp"
#include "DiagUtils.hpp"

namespace tket {

/**
 * Check whether there are any qubits which only requires a single
 * qubit Clifford to make all Paulis I or Z
 */
void check_easy_diagonalise(
    std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> &qubits, Circuit &circ);

/**
 * Given two qubits, attempt to find a basis in which a single CX will
 * make the Paulis on one of qubits fully diagonal
 */
std::optional<std::pair<Pauli, Pauli>> check_pair_compatibility(
    const Qubit &qb1, const Qubit &qb2,
    const std::list<std::pair<QubitPauliTensor, Expr>> &gadgets);

/**
 * Diagonalise a qubit greedily by finding the Pauli Gadget with
 * the lowest residual support over the non-diagonal qubits, and apply
 * single qubit Cliffords and CXs to make it a `ZIII...I` string
 */
void greedy_diagonalise(
    const std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> &qubits, Conjugations &conjugations, Circuit &circ,
    CXConfigType cx_config);

/**
 * Diagonalise a mutually commuting set of Pauli strings. Modifies the
 * list of Pauli strings in place, and returns the Clifford circuit
 * required to generate the initial set.
 */
Circuit mutual_diagonalise(
    std::list<std::pair<QubitPauliTensor, Expr>> &gadgets,
    std::set<Qubit> qubits, CXConfigType cx_config);
/**
 * Applies Clifford conjugations to a QubitPauliTensor
 */
void apply_conjugations(
    QubitPauliTensor &qps, const Conjugations &conjugations);

/**
 *  Given two qubits on which to conjugate a CX gate, try to conjugate with a
 * XXPhase3 instead. If successful, undoes conjugations that must be undone and
 * replaces it with XXPhase3 conjugation. Returns true if successful and false
 * otherwise.
 */
bool conjugate_with_xxphase3(
    const Qubit &qb_a, const Qubit &qb_b, Conjugations &conjugations,
    Circuit &cliff_circ);

}  // namespace tket
